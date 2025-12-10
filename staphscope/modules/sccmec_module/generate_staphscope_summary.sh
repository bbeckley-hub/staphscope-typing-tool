#!/bin/bash

# StaphScope SCCmec Enhanced Summary Generator
# Enhanced to handle K-mer predictions and detailed confidence levels
# Usage: ./generate_staphscope_summary.sh

echo "=================================================="
echo "    STAPHSCOPE Enhanced Summary Report Generator"
echo "=================================================="

# ASCII Art for StaphScope
print_ascii_art() {
    cat << 'EOF'
███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
EOF
}

print_ascii_art

# Output files
TSV_SUMMARY="staphscope_summary.tsv"
HTML_SUMMARY="staphscope_summary.html"
DETAILED_CSV="staphscope_detailed_results.csv"

# Initialize output files with headers
echo -e "Sample_Name\tSCCmec_Type\tMRSA_Status\tConfidence\tPrediction_Method\tKmer_Top_Hits\tEvidence" > "$TSV_SUMMARY"
echo "Sample_Name,SCCmec_Type,MRSA_Status,Confidence,Kmer_Top_Hit_1,Kmer_Coverage_1,Kmer_Top_Hit_2,Kmer_Coverage_2,Evidence,Prediction_Method,Notes" > "$DETAILED_CSV"

# Find all test directories
echo "Scanning for SCCmec results directories..."
s_dirs=$(find . -maxdepth 1 -type d -name "s_*" | sort)

if [ -z "$s_dirs" ]; then
    echo "No s directories found! Please run SCCmec analysis first."
    exit 1
fi

echo "Found $(echo "$s_dirs" | wc -l) result directories"
echo

# Process each test directory
total_processed=0
results_data=()
detailed_results=()

# Statistics counters
mrsa_count=0
mssa_count=0
sccmec_types_count=0
low_confidence_mrsa=0
kmer_based_predictions=0
gene_based_predictions=0
no_type_assigned=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to parse K-mer results from MyKmerFinder output
parse_kmer_results() {
    local kmer_file="$1"
    local top_hits=""
    local coverages=""
    
    if [ -f "$kmer_file" ] && [ -s "$kmer_file" ]; then
        # Get top 2 hits sorted by template coverage (7th column, index 6)
        local hits=$(tail -n +2 "$kmer_file" | head -10 | sort -t$'\t' -k7,7nr | head -2)
        
        local count=0
        while IFS=$'\t' read -r template score expected z pval qcov tcov depth kmers desc; do
            if [ -n "$template" ]; then
                count=$((count + 1))
                # Extract just the type name (before | if exists)
                local type_name=$(echo "$template" | cut -d'|' -f1)
                local coverage=$(echo "$tcov" | sed 's/%//' | xargs)
                
                if [ $count -eq 1 ]; then
                    top_hits="$type_name"
                    coverages="$coverage"
                elif [ $count -eq 2 ]; then
                    top_hits="$top_hits/$type_name"
                    coverages="$coverages/$coverage"
                fi
            fi
        done <<< "$hits"
    fi
    
    echo "$top_hits|$coverages"
}

# Function to parse detailed JSON report
parse_json_report() {
    local json_file="$1"
    local confidence="MEDIUM"
    local evidence=""
    local typing_method="gene_based"
    
    if [ -f "$json_file" ]; then
        # Try to extract confidence from JSON
        confidence=$(jq -r '.mrsa_status.confidence // "MEDIUM"' "$json_file" 2>/dev/null || echo "MEDIUM")
        
        # Try to extract evidence
        evidence=$(jq -r '.mrsa_status.evidence | join("; ")' "$json_file" 2>/dev/null || echo "")
        
        # Try to extract typing method
        typing_method=$(jq -r '.sccmec_typing.typing_method // "gene_based"' "$json_file" 2>/dev/null || echo "gene_based")
        
        # If confidence extraction failed, check if file mentions confidence
        if [[ "$confidence" == "null" || "$confidence" == "" ]]; then
            if grep -q "LOW" "$json_file"; then
                confidence="LOW"
            elif grep -q "HIGH" "$json_file"; then
                confidence="HIGH"
            elif grep -q "VERY_HIGH" "$json_file"; then
                confidence="VERY_HIGH"
            fi
        fi
        
        # Clean up evidence
        evidence=$(echo "$evidence" | sed 's/\"//g' | sed 's/;/ /g')
    fi
    
    echo "$confidence|$evidence|$typing_method"
}

# Function to check if evidence is sufficient for SCCmec assignment
has_sufficient_evidence() {
    local evidence="$1"
    local mrsa_status="$2"
    
    # For MRSA isolates, we need at least mecA detection
    if [[ "$mrsa_status" == "MRSA" ]]; then
        if [[ "$evidence" == *"mecA_detected"* ]] || [[ "$evidence" == *"mecA detected"* ]] || [[ "$evidence" == *"mecC_detected"* ]]; then
            return 0  # true - sufficient evidence
        fi
    fi
    
    # For MSSA, no SCCmec type should be assigned
    return 1  # false - insufficient evidence
}

for dir in $s_dirs; do
    if [ ! -d "$dir" ]; then
        continue
    fi
    
    sample_name="${dir#./s_}"
    results_file="$dir/sccmec_detailed_results.txt"
    json_file="$dir/sccmec_enhanced_report.json"
    kmer_file="$dir/results_MyKmerFinder.txt"
    
    echo -e "${CYAN}Processing: $sample_name${NC}"
    
    if [ ! -f "$results_file" ] && [ ! -f "$json_file" ]; then
        echo -e "  ${YELLOW}⚠️  No results files found, skipping...${NC}"
        continue
    fi
    
    # Initialize variables
    sccmec_type="Not Detected"
    original_sccmec_type=""
    mrsa_status="MSSA"
    confidence="MEDIUM"
    kmer_top_hits=""
    kmer_coverages=""
    evidence=""
    typing_method="gene_based"
    prediction_method=""
    notes=""
    
    # Try to parse JSON report first (more reliable)
    if [ -f "$json_file" ]; then
        echo -e "  📊 Using enhanced JSON report"
        
        # Parse JSON for MRSA status
        if jq -e '.mrsa_status.classification' "$json_file" >/dev/null 2>&1; then
            mrsa_status_raw=$(jq -r '.mrsa_status.classification' "$json_file")
            if [[ "$mrsa_status_raw" == "CONFIRMED_MRSA" ]]; then
                mrsa_status="MRSA"
                mrsa_count=$((mrsa_count + 1))
            else
                mrsa_status="MSSA"
                mssa_count=$((mssa_count + 1))
            fi
        fi
        
        # Parse JSON for SCCmec type
        if jq -e '.sccmec_typing.primary_type' "$json_file" >/dev/null 2>&1; then
            original_sccmec_type=$(jq -r '.sccmec_typing.primary_type' "$json_file")
            if [[ "$original_sccmec_type" != "NA" && "$original_sccmec_type" != "null" ]]; then
                sccmec_type="$original_sccmec_type"
            else
                original_sccmec_type="Not Detected"
            fi
        fi
        
        # Parse JSON for confidence, evidence, and typing method
        json_results=$(parse_json_report "$json_file")
        confidence=$(echo "$json_results" | cut -d'|' -f1)
        evidence=$(echo "$json_results" | cut -d'|' -f2)
        typing_method=$(echo "$json_results" | cut -d'|' -f3)
        
        # Check for k-mer prediction in JSON
        if jq -e '.kmer_prediction' "$json_file" >/dev/null 2>&1; then
            kmer_used=$(jq -r '.kmer_prediction.should_use_kmer' "$json_file")
            if [[ "$kmer_used" == "true" ]]; then
                notes="K-mer based prediction used"
            fi
        fi
        
    else
        # Fall back to parsing text file
        echo -e "  📄 Using text results file"
        
        # Parse the results file
        in_mrsa_section=false
        in_sccmec_section=false
        while IFS= read -r line; do
            line_clean=$(echo "$line" | sed 's/^[ \t]*//;s/[ \t]*$//')
            
            # Check for MRSA/MSSA
            if [[ "$line_clean" == *"MRSA isolate"* ]]; then
                mrsa_status="MRSA"
                mrsa_count=$((mrsa_count + 1))
            elif [[ "$line_clean" == *"MSSA isolate"* ]]; then
                mrsa_status="MSSA"
                mssa_count=$((mssa_count + 1))
            
            # Extract SCCmec type
            elif [[ "$line_clean" == *"Prediction based on genes:"* ]] && [[ "$line_clean" != *"Alert!"* ]]; then
                original_sccmec_type=$(echo "$line_clean" | sed 's/Prediction based on genes: //')
                if [[ "$original_sccmec_type" != "No SCCmec element was detected"* ]]; then
                    sccmec_type="$original_sccmec_type"
                else
                    original_sccmec_type="Not Detected"
                fi
            
            # Check for confidence levels
            elif [[ "$line_clean" == *"Confidence: LOW"* ]]; then
                confidence="LOW"
            elif [[ "$line_clean" == *"Confidence: HIGH"* ]]; then
                confidence="HIGH"
            elif [[ "$line_clean" == *"Confidence: VERY_HIGH"* ]]; then
                confidence="VERY_HIGH"
            
            # Extract evidence
            elif [[ "$line_clean" == *"mecA_detected"* || "$line_clean" == *"mecA detected"* ]]; then
                evidence="$evidence mecA_detected"
            elif [[ "$line_clean" == *"Evidence:"* ]]; then
                evidence_start=true
            elif [[ -n "$evidence_start" && -n "$line_clean" ]]; then
                evidence="$evidence $line_clean"
            fi
            
        done < "$results_file"
        
        # Clean up evidence
        evidence=$(echo "$evidence" | sed 's/^ *//;s/ *$//' | tr -s ' ')
    fi
    
    # Parse K-mer results (for ALL samples)
    if [ -f "$kmer_file" ]; then
        kmer_results=$(parse_kmer_results "$kmer_file")
        if [[ "$kmer_results" != "|" ]]; then
            kmer_top_hits=$(echo "$kmer_results" | cut -d'|' -f1)
            kmer_coverages=$(echo "$kmer_results" | cut -d'|' -f2)
        fi
    fi
    
    # Determine prediction method based on your rules
    if [[ "$confidence" == "HIGH" || "$confidence" == "VERY_HIGH" ]]; then
        prediction_method="Gene-based"
        gene_based_predictions=$((gene_based_predictions + 1))
        
        # For gene-based predictions, use the gene-based type
        if [[ "$sccmec_type" == "Not Detected" || "$sccmec_type" == "NA" ]]; then
            # Should not happen for high confidence, but just in case
            sccmec_type="Gene-based type not found"
        fi
        
    elif [[ "$confidence" == "LOW" || "$confidence" == "MEDIUM" ]]; then
        # Check if we have sufficient evidence
        if has_sufficient_evidence "$evidence" "$mrsa_status"; then
            # Low confidence WITH evidence = K-mer based
            prediction_method="K-mer based"
            kmer_based_predictions=$((kmer_based_predictions + 1))
            
            # Use k-mer top hit as SCCmec type if gene-based didn't find one
            if [[ "$sccmec_type" == "Not Detected" || "$sccmec_type" == "NA" || "$sccmec_type" == "" ]]; then
                if [[ -n "$kmer_top_hits" ]]; then
                    # Take first k-mer hit as the type
                    sccmec_type=$(echo "$kmer_top_hits" | cut -d'/' -f1)
                    notes="K-mer based prediction (low confidence)"
                fi
            fi
            
        else
            # Low confidence WITHOUT evidence = NO SCCmec type assigned
            prediction_method="Insufficient evidence"
            sccmec_type="Not Assigned"
            no_type_assigned=$((no_type_assigned + 1))
            notes="Low confidence without sufficient evidence"
        fi
        
        if [[ "$confidence" == "LOW" && "$mrsa_status" == "MRSA" ]]; then
            low_confidence_mrsa=$((low_confidence_mrsa + 1))
        fi
    fi
    
    # Count SCCmec types (excluding "Not Detected" and "Not Assigned")
    if [[ "$sccmec_type" != "Not Detected" && "$sccmec_type" != "Not Assigned" && "$sccmec_type" != "NA" && "$sccmec_type" != "Gene-based type not found" ]]; then
        sccmec_types_count=$((sccmec_types_count + 1))
    fi
    
    # Format k-mer hits for display (for ALL samples, including gene-based)
    if [[ -n "$kmer_top_hits" ]]; then
        kmer_hit1=$(echo "$kmer_top_hits" | cut -d'/' -f1)
        kmer_cov1=$(echo "$kmer_coverages" | cut -d'/' -f1)
        kmer_hit2=$(echo "$kmer_top_hits" | cut -d'/' -f2 2>/dev/null)
        kmer_cov2=$(echo "$kmer_coverages" | cut -d'/' -f2 2>/dev/null)
        
        if [[ -n "$kmer_hit2" ]]; then
            kmer_display="$kmer_hit1/$kmer_hit2 ($kmer_cov1/$kmer_cov2%)"
        else
            kmer_display="$kmer_hit1 ($kmer_cov1%)"
        fi
    else
        kmer_display="No k-mer hits"
    fi
    
    # Add to TSV
    echo -e "${sample_name}\t${sccmec_type}\t${mrsa_status}\t${confidence}\t${prediction_method}\t${kmer_display}\t${evidence}" >> "$TSV_SUMMARY"
    
    # Add to detailed CSV
    echo "${sample_name},${sccmec_type},${mrsa_status},${confidence},${kmer_hit1},${kmer_cov1},${kmer_hit2},${kmer_cov2},\"${evidence}\",${prediction_method},\"${notes}\"" >> "$DETAILED_CSV"
    
    # Store for HTML
    results_data+=("$sample_name" "$sccmec_type" "$mrsa_status" "$confidence" "$prediction_method" "$kmer_display" "$evidence")
    
    # Color-coded output
    status_color="$GREEN"
    [[ "$mrsa_status" == "MRSA" ]] && status_color="$RED"
    
    confidence_color="$YELLOW"
    [[ "$confidence" == "HIGH" ]] && confidence_color="$GREEN"
    [[ "$confidence" == "VERY_HIGH" ]] && confidence_color="$BLUE"
    [[ "$confidence" == "LOW" ]] && confidence_color="$RED"
    
    method_color="$CYAN"
    [[ "$prediction_method" == "K-mer based" ]] && method_color="$MAGENTA"
    [[ "$prediction_method" == "Insufficient evidence" ]] && method_color="$YELLOW"
    
    echo -e "  ${status_color}✅ $mrsa_status${NC} | ${confidence_color}$confidence${NC} | ${method_color}$prediction_method${NC}"
    echo -e "  📋 Type: $sccmec_type"
    if [[ -n "$kmer_display" && "$kmer_display" != "No k-mer hits" ]]; then
        echo -e "    🔍 K-mer hits: $kmer_display"
    fi
    if [[ -n "$evidence" ]]; then
        echo -e "    📝 Evidence: $evidence"
    fi
    
    total_processed=$((total_processed + 1))
done

# Generate HTML Summary Report
generate_html_report() {
    cat << EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - Enhanced SCCmec Summary Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1800px;
            margin: 0 auto;
        }
        
        .header {
            text-align: center;
            margin-bottom: 30px;
        }
        
        .ascii-container {
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 255, 0, 0.3);
        }
        
        .ascii-art {
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.1;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 10px rgba(0, 255, 0, 0.5);
            text-align: center;
            overflow-x: auto;
        }
        
        .quote-container {
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }
        
        .quote-text {
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
            color: #ffffff;
        }
        
        .quote-author {
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }
        
        .summary-section {
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }
        
        .summary-section h2 {
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }
        
        .results-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 12px;
        }
        
        .results-table th {
            background: #1e40af;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
        }
        
        .results-table td {
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }
        
        .results-table tr:hover {
            background: #f3f4f6;
        }
        
        .status-mrsa {
            background: #dc2626;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 11px;
        }
        
        .status-mssa {
            background: #16a34a;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 11px;
        }
        
        .confidence-high {
            background: #16a34a;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .confidence-medium {
            background: #f59e0b;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .confidence-low {
            background: #dc2626;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .confidence-very-high {
            background: #1e40af;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .method-gene {
            background: #10b981;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .method-kmer {
            background: #8b5cf6;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .method-insufficient {
            background: #6b7280;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 10px;
        }
        
        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }
        
        .metric-card {
            background: linear-gradient(135deg, #8b5cf6 0%, #6d28d9 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
            text-align: center;
        }
        
        .metric-label {
            font-size: 12px;
            opacity: 0.9;
            margin-bottom: 5px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        .metric-value {
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }
        
        .metric-subtitle {
            font-size: 10px;
            opacity: 0.8;
        }
        
        .footer {
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }
        
        .timestamp {
            color: #fbbf24;
            font-weight: bold;
        }
        
        .authorship {
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }
        
        .evidence-cell {
            font-size: 10px;
            color: #4b5563;
            max-width: 300px;
            overflow: visible;
            white-space: normal;
            word-wrap: break-word;
        }
        
        .no-assignment {
            color: #dc2626;
            font-weight: bold;
            font-style: italic;
        }
        
        .kmer-cell {
            font-size: 10px;
            color: #374151;
        }
        
        @media (max-width: 768px) {
            .results-table {
                font-size: 10px;
            }
            .ascii-art {
                font-size: 8px;
            }
            .metrics-grid {
                grid-template-columns: repeat(2, 1fr);
            }
            .evidence-cell {
                max-width: 150px;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
                </div>
            </div>
            
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText"></div>
                <div class="quote-author" id="quoteAuthor"></div>
            </div>
        </div>
        
        <div class="summary-section">
            <h2>📊 SCCmec Analysis Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Total Samples</div>
                    <div class="metric-value">$total_processed</div>
                    <div class="metric-subtitle">Processed Successfully</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MRSA Isolates</div>
                    <div class="metric-value">$mrsa_count</div>
                    <div class="metric-subtitle">
                        $low_confidence_mrsa low confidence
                    </div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">MSSA Isolates</div>
                    <div class="metric-value">$mssa_count</div>
                    <div class="metric-subtitle">No mecA/mecC detected</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">SCCmec Types</div>
                    <div class="metric-value">$sccmec_types_count</div>
                    <div class="metric-subtitle">
                        Gene-based: $gene_based_predictions
                    </div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">K-mer Based</div>
                    <div class="metric-value">$kmer_based_predictions</div>
                    <div class="metric-subtitle">Low confidence with evidence</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Not Assigned</div>
                    <div class="metric-value">$no_type_assigned</div>
                    <div class="metric-subtitle">Low confidence, no evidence</div>
                </div>
            </div>
        </div>
        
        <div class="summary-section">
            <h2>🧬 Detailed Results</h2>
            <div style="max-height: 800px; overflow-y: auto;">
                <table class="results-table">
                    <thead>
                        <tr>
                            <th>Sample Name</th>
                            <th>SCCmec Type</th>
                            <th>MRSA Status</th>
                            <th>Confidence</th>
                            <th>Prediction Method</th>
                            <th>K-mer Top Hits</th>
                            <th>Evidence</th>
                        </tr>
                    </thead>
                    <tbody>
EOF

    # Add table rows from results data
    for ((i=0; i<${#results_data[@]}; i+=7)); do
        sample="${results_data[i]}"
        sccmec="${results_data[i+1]}"
        status="${results_data[i+2]}"
        confidence="${results_data[i+3]}"
        method="${results_data[i+4]}"
        kmer_hits="${results_data[i+5]}"
        evidence="${results_data[i+6]}"
        
        # Determine status class
        status_class="status-mrsa"
        [[ "$status" == "MSSA" ]] && status_class="status-mssa"
        
        # Determine confidence class
        confidence_class="confidence-medium"
        [[ "$confidence" == "HIGH" ]] && confidence_class="confidence-high"
        [[ "$confidence" == "VERY_HIGH" ]] && confidence_class="confidence-very-high"
        [[ "$confidence" == "LOW" ]] && confidence_class="confidence-low"
        
        # Determine method class
        method_class="method-gene"
        [[ "$method" == "K-mer based" ]] && method_class="method-kmer"
        [[ "$method" == "Insufficient evidence" ]] && method_class="method-insufficient"
        
        # Format SCCmec type cell
        if [[ "$sccmec" == "Not Assigned" ]]; then
            sccmec_cell="<span class='no-assignment'>$sccmec</span>"
        else
            sccmec_cell="$sccmec"
        fi
        
        echo "<tr>"
        echo "<td><strong>$sample</strong></td>"
        echo "<td>$sccmec_cell</td>"
        echo "<td><span class=\"$status_class\">$status</span></td>"
        echo "<td><span class=\"$confidence_class\">$confidence</span></td>"
        echo "<td><span class=\"$method_class\">$method</span></td>"
        echo "<td class=\"kmer-cell\">$kmer_hits</td>"
        echo "<td class=\"evidence-cell\">$evidence</td>"
        echo "</tr>"
    done

    cat << EOF
                    </tbody>
                </table>
            </div>
            <div style="margin-top: 20px; padding: 10px; background: #f8fafc; border-radius: 6px; font-size: 11px; color: #4b5563;">
                <strong>Interpretation Rules:</strong><br>
                1. <span class="method-gene" style="margin-left: 5px;">Gene-based</span> = High/Very High confidence predictions<br>
                2. <span class="method-kmer" style="margin-left: 5px;">K-mer based</span> = Low confidence WITH evidence (mecA/mecC detected)<br>
                3. <span class="method-insufficient" style="margin-left: 5px;">Insufficient evidence</span> = Low confidence WITHOUT evidence → No SCCmec type assigned<br>
                4. <span class="no-assignment" style="margin-left: 5px;">Not Assigned</span> = Low confidence, no evidence (shown in red)<br>
                5. All samples show k-mer hits for reference (gene-based and k-mer based)
            </div>
        </div>
        
        <div class="footer">
            <p><strong>STAPHSCOPE</strong> - Enhanced SCCmec Cassette Analysis Summary</p>
            <p class="timestamp">Generated: $(date "+%Y-%m-%d %H:%M:%S")</p>
            <p style="margin-top: 10px; font-size: 12px;">
                Total samples: $total_processed | 
                Gene-based: $gene_based_predictions | 
                K-mer based: $kmer_based_predictions |
                Not assigned: $no_type_assigned
            </p>
            
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley</p>
                <p>GitHub: <a href="https://github.com/bbeckley-hub" style="color: #fbbf24;">bbeckley-hub</a></p>
                <p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color: #fbbf24;">brownbeckley94@gmail.com</a></p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>

    <script>
        const quotes = [
            { text: "The important thing is not to stop questioning. Curiosity has its own reason for existing.", author: "Albert Einstein" },
            { text: "Science is not only a disciple of reason but also one of romance and passion.", author: "Stephen Hawking" },
            { text: "Somewhere, something incredible is waiting to be known.", author: "Carl Sagan" },
            { text: "The good thing about science is that it's true whether or not you believe in it.", author: "Neil deGrasse Tyson" },
            { text: "In science, there are no shortcuts to truth.", author: "Karl Popper" },
            { text: "Science knows no country, because knowledge belongs to humanity.", author: "Louis Pasteur" },
            { text: "The science of today is the technology of tomorrow.", author: "Edward Teller" },
            { text: "Nothing in life is to be feared, it is only to be understood.", author: "Marie Curie" },
            { text: "Research is what I'm doing when I don't know what I'm doing.", author: "Wernher von Braun" },
            { text: "The universe is not required to be in perfect harmony with human ambition.", author: "Carl Sagan" }
        ];

        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');

        function getRandomQuote() {
            return quotes[Math.floor(Math.random() * quotes.length)];
        }

        function displayQuote() {
            quoteContainer.style.opacity = '0';
            
            setTimeout(() => {
                const quote = getRandomQuote();
                quoteText.textContent = "\"" + quote.text + "\"";
                quoteAuthor.textContent = "— " + quote.author;
                quoteContainer.style.opacity = '1';
            }, 500);
        }

        displayQuote();
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>
EOF
}

# Generate HTML report
generate_html_report > "$HTML_SUMMARY"

echo
echo "=================================================="
echo "           Enhanced Summary Generation Completed!"
echo "=================================================="
echo -e "${GREEN}📁 TSV Summary: $TSV_SUMMARY${NC}"
echo -e "${GREEN}📊 Detailed CSV: $DETAILED_CSV${NC}"
echo -e "${GREEN}🌐 HTML Summary: $HTML_SUMMARY${NC}"
echo -e "${CYAN}📈 Total samples processed: $total_processed${NC}"
echo
echo -e "${BLUE}Summary Statistics:${NC}"
echo "-------------------"
echo -e "${RED}MRSA isolates: $mrsa_count${NC}"
echo -e "  └ Low confidence MRSA: $low_confidence_mrsa"
echo -e "${GREEN}MSSA isolates: $mssa_count${NC}"
echo -e "${YELLOW}SCCmec types assigned: $sccmec_types_count${NC}"
echo -e "${MAGENTA}Gene-based predictions: $gene_based_predictions${NC}"
echo -e "${MAGENTA}K-mer based predictions: $kmer_based_predictions${NC}"
echo -e "${RED}Not assigned (no evidence): $no_type_assigned${NC}"
echo
echo -e "${CYAN}✅ Analysis completed successfully!${NC}"
echo

# Display a sample of the summary table
if [ $total_processed -gt 0 ]; then
    echo -e "${BLUE}Sample of Results:${NC}"
    echo "---------------------------------------------------------------------------------------------------------------"
    head -5 "$TSV_SUMMARY" | column -t -s $'\t'
    echo "---------------------------------------------------------------------------------------------------------------"
fi

echo
print_ascii_art
