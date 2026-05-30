from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
import subprocess
import json
import os

app = FastAPI()

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class AnalysisRequest(BaseModel):
    dna: str
    pattern: str

@app.get("/")
async def get_dashboard():
    # Serves the dashboard.html at http://127.0.0.1:8001/
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "dashboard.html")
    return FileResponse(file_path)

@app.post("/analyze")
async def analyze_dna(request: AnalysisRequest):
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        executable = os.path.join(base_dir, "GenomeToJSON.exe")
        if not os.path.exists(executable):
            raise HTTPException(status_code=500, detail="C++ Engine not found. Please compile it first.")

        result = subprocess.run(
            [executable, request.dna, request.pattern],
            capture_output=True,
            text=True,
            check=True
        )
        
        output = result.stdout.strip()
        
        try:
            start_idx = output.find('{')
            end_idx = output.rfind('}') + 1
            if start_idx != -1 and end_idx != 0:
                json_str = output[start_idx:end_idx]
                data = json.loads(json_str)
                return data
            else:
                raise ValueError("No JSON found in output")
        except Exception as e:
            raise HTTPException(status_code=500, detail="Failed to parse engine output")

    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=500, detail=f"Engine error: {e.stderr}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    print("\nDNA Analysis Server starting...")
    print("Open your browser at: http://127.0.0.1:8001")
    uvicorn.run(app, host="127.0.0.1", port=8001)
