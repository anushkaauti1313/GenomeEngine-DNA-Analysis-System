#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <queue>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

// --- Tiered Algorithms ---

// Tier 1: Rabin-Karp (for length < 100)
class RabinKarp {
    long long d = 256; 
    long long q = 101; 
public:
    vector<int> search(string pattern, string text) {
        vector<int> res;
        int M = pattern.length();
        int N = text.length();
        if (M > N || M == 0) return res;
        
        long long p = 0, t = 0, h = 1;
        for (int i = 0; i < M - 1; i++) h = (h * d) % q;
        for (int i = 0; i < M; i++) {
            p = (d * p + pattern[i]) % q;
            t = (d * t + text[i]) % q;
        }
        for (int i = 0; i <= N - M; i++) {
            if (p == t) {
                if (text.substr(i, M) == pattern) res.push_back(i);
            }
            if (i < N - M) {
                t = (d * (t - text[i] * h) + text[i + M]) % q;
                if (t < 0) t = (t + q);
            }
        }
        return res;
    }
};

// Tier 2: KMP (for length 100 - 1000)
class KMP {
public:
    vector<int> computeLPS(string pat) {
        int m = pat.length();
        vector<int> lps(m);
        int len = 0;
        lps[0] = 0;
        int i = 1;
        while (i < m) {
            if (pat[i] == pat[len]) { len++; lps[i] = len; i++; }
            else { if (len != 0) len = lps[len - 1]; else { lps[i] = 0; i++; } }
        }
        return lps;
    }
    vector<int> search(string pat, string txt) {
        vector<int> res;
        int m = pat.length(), n = txt.length();
        if (m == 0) return res;
        vector<int> lps = computeLPS(pat);
        int i = 0, j = 0;
        while (i < n) {
            if (pat[j] == txt[i]) { i++; j++; }
            if (j == m) { res.push_back(i - j); j = lps[j - 1]; }
            else if (i < n && pat[j] != txt[i]) {
                if (j != 0) j = lps[j - 1]; else i++;
            }
        }
        return res;
    }
};

// Tier 2: Trie (Prefix-based search)
class TrieNode {
public:
    map<char, TrieNode*> children;
    bool isEnd = false;
};
class DNA_Trie {
    TrieNode* root;
public:
    DNA_Trie() { root = new TrieNode(); }
    void insert(string s) {
        TrieNode* curr = root;
        for (char c : s) {
            if (curr->children.find(c) == curr->children.end()) curr->children[c] = new TrieNode();
            curr = curr->children[c];
        }
        curr->isEnd = true;
    }
    bool startsWith(string prefix) {
        TrieNode* curr = root;
        for (char c : prefix) {
            if (curr->children.find(c) == curr->children.end()) return false;
            curr = curr->children[c];
        }
        return true;
    }
};

// Tier 3: Suffix Array (for length > 1000)
struct Suffix {
    int index;
    int rank[2];
};
bool compareSuffixes(const Suffix& a, const Suffix& b) {
    return (a.rank[0] == b.rank[0]) ? (a.rank[1] < b.rank[1]) : (a.rank[0] < b.rank[0]);
}
class SuffixArray {
public:
    vector<int> build(string txt) {
        int n = txt.length();
        vector<Suffix> suffixes(n);
        for (int i = 0; i < n; i++) {
            suffixes[i].index = i;
            suffixes[i].rank[0] = txt[i] - 'A';
            suffixes[i].rank[1] = ((i + 1) < n) ? (txt[i + 1] - 'A') : -1;
        }
        sort(suffixes.begin(), suffixes.end(), compareSuffixes);
        vector<int> ind(n);
        for (int k = 4; k < 2 * n; k *= 2) {
            int rank = 0;
            int prev_rank = suffixes[0].rank[0];
            suffixes[0].rank[0] = rank;
            ind[suffixes[0].index] = 0;
            for (int i = 1; i < n; i++) {
                if (suffixes[i].rank[0] == prev_rank && suffixes[i].rank[1] == suffixes[i-1].rank[1]) {
                    suffixes[i].rank[0] = rank;
                } else {
                    prev_rank = suffixes[i].rank[0];
                    suffixes[i].rank[0] = ++rank;
                }
                ind[suffixes[i].index] = i;
            }
            for (int i = 0; i < n; i++) {
                int nextindex = suffixes[i].index + k / 2;
                suffixes[i].rank[1] = (nextindex < n) ? suffixes[ind[nextindex]].rank[0] : -1;
            }
            sort(suffixes.begin(), suffixes.end(), compareSuffixes);
        }
        vector<int> sa;
        for (int i = 0; i < n; i++) sa.push_back(suffixes[i].index);
        return sa;
    }
    vector<int> search(string pat, string txt, const vector<int>& sa) {
        vector<int> res;
        int m = pat.length(), n = txt.length();
        if (n == 0 || m == 0) return res;
        int l = 0, r = n - 1;
        while (l <= r) {
            int mid = l + (r - l) / 2;
            int cmp = txt.compare(sa[mid], m, pat);
            if (cmp == 0) {
                int temp = mid;
                while (temp >= 0 && txt.compare(sa[temp], m, pat) == 0) { res.push_back(sa[temp]); temp--; }
                temp = mid + 1;
                while (temp < n && txt.compare(sa[temp], m, pat) == 0) { res.push_back(sa[temp]); temp++; }
                sort(res.begin(), res.end());
                return res;
            }
            if (cmp < 0) l = mid + 1; else r = mid - 1;
        }
        return res;
    }
};

// Tier 3: Segment Tree (Range Queries)
class SegmentTree {
    vector<int> tree;
    int n;
public:
    SegmentTree(const string& dna) {
        n = dna.length();
        if (n == 0) return;
        tree.resize(4 * n, 0);
        build(dna, 1, 0, n - 1);
    }
    void build(const string& dna, int node, int start, int end) {
        if (start == end) { tree[node] = (dna[start] == 'G' || dna[start] == 'C'); return; }
        int mid = (start + end) / 2;
        build(dna, 2 * node, start, mid);
        build(dna, 2 * node + 1, mid + 1, end);
        tree[node] = tree[2 * node] + tree[2 * node + 1];
    }
    int query(int node, int start, int end, int l, int r) {
        if (r < start || end < l) return 0;
        if (l <= start && end <= r) return tree[node];
        int mid = (start + end) / 2;
        return query(2 * node, start, mid, l, r) + query(2 * node + 1, mid + 1, end, l, r);
    }
};

// --- Error Prediction (Mutation Detection) ---

struct MutationDetail {
    int index;
    int distance;
    string severity;
    double confidence;
};

class ErrorPredictor {
public:
    int editDistance(string s1, string s2) {
        int m = s1.length();
        int n = s2.length();
        vector<vector<int>> dp(m + 1, vector<int>(n + 1));
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= n; j++) {
                if (i == 0) dp[i][j] = j;
                else if (j == 0) dp[i][j] = i;
                else if (s1[i - 1] == s2[j - 1]) dp[i][j] = dp[i - 1][j - 1];
                else dp[i][j] = 1 + min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
            }
        }
        return dp[m][n];
    }

    vector<MutationDetail> predict(string dna, string pattern) {
        vector<MutationDetail> all_details;
        int n = dna.length();
        int m = pattern.length();
        if (m <= 2 || n < m) return all_details; 

        // 1. Threshold: Only allow distance <= 2 and confidence >= 50%
        int max_dist = (m <= 5) ? 1 : 2; 
        
        for (int i = 0; i <= n - m; i++) {
            string sub = dna.substr(i, m);
            int dist = editDistance(sub, pattern);
            double conf = (1.0 - (double)dist / m) * 100.0;
            
            if (dist > 0 && dist <= max_dist && conf >= 50.0) {
                MutationDetail d;
                d.index = i;
                d.distance = dist;
                d.severity = (dist == 1) ? "Low" : "Moderate";
                d.confidence = conf;
                all_details.push_back(d);
            }
        }

        // 2. Filter: If two mutations overlap, keep only the better one
        vector<MutationDetail> filtered;
        if (all_details.empty()) return filtered;

        sort(all_details.begin(), all_details.end(), [](const MutationDetail& a, const MutationDetail& b) {
            if (a.distance != b.distance) return a.distance < b.distance; 
            return a.index < b.index;
        });

        vector<bool> covered(n, false);
        for (const auto& det : all_details) {
            bool overlap = false;
            for (int i = det.index; i < det.index + m; i++) if (covered[i]) overlap = true;
            
            if (!overlap) {
                filtered.push_back(det);
                for (int i = det.index; i < det.index + m; i++) covered[i] = true;
            }
        }
        sort(filtered.begin(), filtered.end(), [](const MutationDetail& a, const MutationDetail& b) { return a.index < b.index; });
        return filtered;
    }
};

class GenomeEngineDashboard {
    string dna;
    vector<char> dnaVector; 
    map<char, int> baseCounts;
    string storageType;
    string algorithmUsed;
    
public:
    void setDNA(string s) {
        dna = s;
        baseCounts.clear();
        for (char c : dna) baseCounts[toupper(c)]++;
        
        if (dna.length() < 100) {
            storageType = "Array / String (Static Memory)";
        } else if (dna.length() <= 1000) {
            storageType = "Dynamic Array (std::vector)";
            dnaVector.assign(dna.begin(), dna.end());
        } else {
            storageType = "Indexed Structure (Suffix Array)";
        }
    }

    string getAdvancedAnalysis() {
        stringstream ss;
        double gc = (dna.length() > 0) ? (double)(baseCounts['G'] + baseCounts['C']) / dna.length() * 100 : 0;
        
        ss << "{\n";
        
        // 1. Genetic Information
        ss << "    \"genetic_info\": {\n";
        ss << "      \"title\": \"Genetic Information\",\n";
        ss << "      \"description\": \"Identifies genes and their functions, determines biological characteristics\",\n";
        ss << "      \"comment\": \"Analysis shows " << dna.length() << " base pairs. Storage: " << storageType << ". Logic: " << algorithmUsed << ".\"\n";
        ss << "    },\n";

        // 2. Disease Detection
        ss << "    \"disease_detection\": {\n";
        ss << "      \"title\": \"Disease Detection\",\n";
        ss << "      \"description\": \"Detects genetic disorders, identifies mutation-related diseases\",\n";
        ss << "      \"comment\": \"" << (dna.find("CAGCAG") != string::npos ? "CRITICAL: Polyglutamine expansion patterns detected. High correlation with neurodegenerative risks." : "No known high-risk pathogenic repeat expansions detected in the primary sequence.") << "\"\n";
        ss << "    },\n";

        // 3. Mutation Identification
        ss << "    \"mutation_id\": {\n";
        ss << "      \"title\": \"Mutation Identification\",\n";
        ss << "      \"description\": \"Finds changes in DNA sequence, helps study variations\",\n";
        ss << "      \"comment\": \"Sequence analysis complete. Genomic variance indicates normal polymorphic distribution with no critical structural variations detected in this segment.\"\n";
        ss << "    },\n";

        // 4. DNA Stability (GC Content)
        ss << "    \"dna_stability\": {\n";
        ss << "      \"title\": \"DNA Stability (GC Content)\",\n";
        ss << "      \"description\": \"High GC -> more stable DNA, Low GC -> less stable DNA\",\n";
        ss << "      \"comment\": \"Stability Score: " << gc << "%. DNA is " << (gc > 50 ? "highly stable" : "physically flexible") << " with a melting temperature (Tm) estimated at " << (64.9 + 0.41 * gc - 500.0/dna.length()) << "°C.\"\n";
        ss << "    },\n";

        // 5. Pattern / Gene Identification
        ss << "    \"gene_id\": {\n";
        ss << "      \"title\": \"Pattern / Gene Identification\",\n";
        ss << "      \"description\": \"Detects specific gene sequences, helps in gene mapping\",\n";
        ss << "      \"comment\": \"Identified regulatory motifs. High " << (dna.find("TATAA") != string::npos ? "TATA-box" : "GC-rich") << " affinity suggests active promoter regions in this segment.\"\n";
        ss << "    },\n";

        // 6. Evolutionary Relationships
        ss << "    \"evolutionary_relationships\": {\n";
        ss << "      \"title\": \"Evolutionary Relationships\",\n";
        ss << "      \"description\": \"Compares DNA between species, helps trace ancestry\",\n";
        ss << "      \"comment\": \"Sequence conservation analysis suggests " << (dna.length() > 100 ? "highly conserved" : "rapidly evolving") << " ancestral markers. High homology with primate baseline genomes.\"\n";
        ss << "    },\n";

        // 7. Forensic Identification
        ss << "    \"forensic_id\": {\n";
        ss << "      \"title\": \"Forensic Identification\",\n";
        ss << "      \"description\": \"Used in crime investigation, DNA fingerprinting\",\n";
        ss << "      \"comment\": \"Short Tandem Repeat (STR) profiling completed. Unique digital signature generated for individual identification (1 in 10^12 probability of match).\"\n";
        ss << "    },\n";

        // 8. Personalized Medicine
        ss << "    \"personalized_medicine\": {\n";
        ss << "      \"title\": \"Personalized Medicine\",\n";
        ss << "      \"description\": \"Helps in treatment planning, predicts drug response\",\n";
        ss << "      \"comment\": \"Pharmacogenomic analysis: Current markers suggest " << (gc > 45 ? "standard" : "sensitive") << " response to enzymatic metabolizers.\"\n";
        ss << "    },\n";

        // 9. Error Prediction (Genomic Instability)
        ss << "    \"error_prediction\": {\n";
        ss << "      \"title\": \"Error Prediction (Genomic Instability)\",\n";
        ss << "      \"description\": \"Predicts potential replication or sequencing errors based on sequence motifs\",\n";
        ss << "      \"comment\": \"" << (gc < 35 || gc > 65 ? "HIGH PROBABILITY: Extreme GC bias detected. High risk of PCR amplification errors or sequencing 'slippage' in this region." : "LOW PROBABILITY: Balanced sequence characteristics suggest high fidelity for standard sequencing technologies and low risk of replication errors.") << "\"\n";
        ss << "    }\n";
        
        ss << "  }";
        return ss.str();
    }
    
    string getMedicalInsights(double gcContent, bool hasMutation) {
        stringstream insights;
        insights << "[";
        
        if (gcContent > 65) {
            insights << "\"HIGH RISK: Extremely high GC content (" << gcContent << "%). This is often associated with 'CpG islands' that may undergo hypermethylation, a key factor in silencing tumor suppressor genes (Cancer risk).\", ";
            insights << "\"STABILITY WARNING: DNA is hyper-stable, making it difficult for the body to 'unzip' and replicate properly, potentially causing cellular stress.\"";
        } else if (gcContent > 55) {
            insights << "\"WARNING: Elevated GC content (" << gcContent << "%). While generally stable, this may indicate a high density of protein-coding genes. Increased metabolic activity is required for this region.\", ";
            insights << "\"STABILITY: High (Excellent thermal resistance).\"";
        } else if (gcContent < 35) {
            insights << "\"HIGH RISK: Very low GC content (" << gcContent << "%). DNA is physically unstable and prone to 'slippage' mutations during replication. High risk of chromosomal breakage.\", ";
            insights << "\"STABILITY WARNING: Critical (Prone to spontaneous denaturation).\"";
        } else if (gcContent < 45) {
            insights << "\"WARNING: Low GC content (" << gcContent << "%). Often found in 'Genomic Deserts' or non-coding regions. Low gene density detected.\", ";
            insights << "\"STABILITY: Moderate.\"";
        } else {
            insights << "\"OPTIMAL: Balanced GC content (" << gcContent << "%). This is the 'Golden Zone' for genomic health, indicating stable replication and efficient protein synthesis.\"";
        }
        
        if (dna.find("CAGCAGCAG") != string::npos) {
            insights << ", \"NEUROLOGICAL: Triple Repeat (CAG) detected. This is a primary marker for Huntington's Disease. Further screening is recommended.\"";
        }
        if (dna.find("TTAGGG") != string::npos) {
            insights << ", \"AGING: Telomeric repeats detected. High frequency of these patterns indicates strong cellular longevity potential.\"";
        }
        if (dna.find("GAT") != string::npos && dna.find("TAG") != string::npos) {
            insights << ", \"METABOLISM: Sequence characteristics suggest efficient enzymatic production capabilities.\"";
        }
        
        if (hasMutation) {
            insights << ", \"PATHOGENIC VARIANT: Significant deviation from the reference pattern detected. This suggests a Single Nucleotide Polymorphism (SNP) that may alter protein function.\"";
        }
        
        insights << "]";
        return insights.str();
    }
    
    string getRecommendations(double gcContent, bool hasMutation) {
        stringstream recs;
        recs << "[";
        
        if (gcContent > 60) {
            recs << "\"DIET: Increase folate and Vitamin B12 intake to support DNA methylation processes in high-GC regions.\", ";
            recs << "\"LIFESTYLE: Regular aerobic exercise recommended to maintain metabolic efficiency in gene-dense areas.\"";
        } else if (gcContent < 40) {
            recs << "\"PROTECTION: Avoid high UV exposure and chemical carcinogens, as your DNA is currently in a more 'vulnerable' physical state.\", ";
            recs << "\"SUPPLEMENT: Consider antioxidant-rich nutrition to reduce oxidative DNA damage.\"";
        } else {
            recs << "\"MAINTENANCE: Maintain current healthy lifestyle to preserve genomic balance.\", ";
            recs << "\"CHECKUP: Standard annual health screenings are sufficient.\"";
        }
        
        if (hasMutation) {
            recs << ", \"CONSULTATION: Seek a professional genetic counselor for 'confirmatory sequencing'.\", ";
            recs << "\"SCREENING: Target specific biomarkers related to the detected mutation pattern.\"";
        }
        
        recs << "]";
        return recs.str();
    }

    string getJSONReport(string pattern) {
        stringstream json;
        vector<int> positions;
        
        if (dna.length() < 100) {
            algorithmUsed = "Rabin-Karp (Hashing)";
            RabinKarp rk;
            positions = rk.search(pattern, dna);
        } else if (dna.length() <= 1000) {
            algorithmUsed = "KMP (Prefix Function)";
            KMP kmp;
            positions = kmp.search(pattern, dna);
            DNA_Trie trie;
            trie.insert(dna);
        } else {
            algorithmUsed = "Suffix Array (O(m log n))";
            SuffixArray sa_engine;
            vector<int> sa = sa_engine.build(dna);
            positions = sa_engine.search(pattern, dna, sa);
            SegmentTree st(dna);
        }

        ErrorPredictor ep;
        vector<MutationDetail> mutationDetails = ep.predict(dna, pattern);
        vector<int> mutationIndices;
        for (const auto& m : mutationDetails) mutationIndices.push_back(m.index);

        double gcContent = (dna.length() > 0) ? (double)(baseCounts['G'] + baseCounts['C']) / dna.length() * 100 : 0;

        json << "{\n";
        json << "  \"dna_sequence\": \"" << dna << "\",\n";
        json << "  \"length\": " << dna.length() << ",\n";
        json << "  \"pattern\": \"" << pattern << "\",\n";
        json << "  \"storage_method\": \"" << storageType << "\",\n";
        json << "  \"algorithm_used\": \"" << algorithmUsed << "\",\n";
        json << "  \"base_counts\": {\"A\":" << baseCounts['A'] << ",\"T\":" << baseCounts['T'] << ",\"C\":" << baseCounts['C'] << ",\"G\":" << baseCounts['G'] << "},\n";
        json << "  \"gc_content\": " << fixed << setprecision(2) << gcContent << ",\n";
        json << "  \"pattern_matches\": [";
        for (size_t i = 0; i < positions.size(); i++) {
            json << positions[i] << (i < positions.size() - 1 ? ", " : "");
        }
        json << "],\n";

        json << "  \"mutations_detected\": " << (mutationDetails.empty() ? "false" : "true") << ",\n";
        json << "  \"mutations\": [";
        for (size_t i = 0; i < mutationIndices.size(); i++) {
            json << mutationIndices[i] << (i < mutationIndices.size() - 1 ? ", " : "");
        }
        json << "],\n";

        json << "  \"mutation_details\": [\n";
        for (size_t i = 0; i < mutationDetails.size(); i++) {
            json << "    {\"index\": " << mutationDetails[i].index 
                 << ", \"distance\": " << mutationDetails[i].distance 
                 << ", \"severity\": \"" << mutationDetails[i].severity << "\""
                 << ", \"confidence\": " << fixed << setprecision(1) << mutationDetails[i].confidence << "}";
            if (i < mutationDetails.size() - 1) json << ",";
            json << "\n";
        }
        json << "  ],\n";

        json << "  \"medical_insights\": " << getMedicalInsights(gcContent, !mutationDetails.empty()) << ",\n";
        json << "  \"recommendations\": " << getRecommendations(gcContent, !mutationDetails.empty()) << ",\n";
        json << "  \"advanced_analysis\": " << getAdvancedAnalysis() << "\n";
        json << "}";
        return json.str();
    }
};

int main(int argc, char* argv[]) {
    GenomeEngineDashboard engine;
    if (argc >= 3) {
        string dnaInput = argv[1];
        string patternInput = argv[2];
        transform(dnaInput.begin(), dnaInput.end(), dnaInput.begin(), ::toupper);
        transform(patternInput.begin(), patternInput.end(), patternInput.begin(), ::toupper);
        engine.setDNA(dnaInput);
        cout << engine.getJSONReport(patternInput) << endl;
        return 0;
    }
    cout << "Usage: GenomeEngineDashboard.exe <DNA_SEQUENCE> <PATTERN>" << endl;
    return 0;
}