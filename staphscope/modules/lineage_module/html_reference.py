# staphscope/modules/lineage_module/html_reference.py
#!/usr/bin/env python3
"""
StaphScope Lineage Reference HTML Generator
Creates a comprehensive reference document with ALL database content
Author: Brown Beckley <brownbeckley94@gmail.com>
"""

import os
from datetime import datetime
from lineage_database import LINEAGE_DATABASE, SPECIALIZED_LINEAGES, GEOGRAPHIC_LINEAGES

def generate_comprehensive_lineage_html(output_path="staphscope_lineage_reference.html"):
    """Generate comprehensive HTML reference with ALL database content"""
    
    lineages = LINEAGE_DATABASE
    specialized = SPECIALIZED_LINEAGES
    geographic = GEOGRAPHIC_LINEAGES
    
    # Calculate statistics
    total_lineages = len(lineages)
    stats = {
        "healthcare": sum(1 for info in lineages.values() if "Healthcare" in info.get("type", "")),
        "community": sum(1 for info in lineages.values() if "Community" in info.get("type", "")),
        "livestock": sum(1 for info in lineages.values() if "Livestock" in info.get("type", "")),
        "high_risk": sum(1 for info in lineages.values() if info.get("risk_level") == "HIGH"),
        "medium_risk": sum(1 for info in lineages.values() if info.get("risk_level") == "MEDIUM"),
        "low_risk": sum(1 for info in lineages.values() if info.get("risk_level") == "LOW"),
        "pvl_positive": len(specialized.get("PVL_POSITIVE", {})),
        "tst_positive": len(specialized.get("TST_POSITIVE", {})),
        "mec_c": len(specialized.get("MEC_C_POSITIVE", {})),
    }
    
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>STAPHSCOPE - Comprehensive S. aureus Lineage Reference</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary: #1e3a8a;
            --secondary: #3b82f6;
            --accent: #f59e0b;
            --danger: #dc2626;
            --success: #16a34a;
            --warning: #f59e0b;
            --dark: #1f2937;
            --light: #f8fafc;
        }}
        
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 50%, #7e22ce 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: var(--light);
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        /* Header Styles */
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0, 0, 0, 0.8);
            padding: 25px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.6);
            border: 3px solid rgba(0, 255, 0, 0.4);
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 14px;
            line-height: 1.3;
            white-space: pre;
            color: #00ff00;
            text-shadow: 0 0 15px rgba(0, 255, 0, 0.7);
            text-align: center;
            letter-spacing: 1px;
        }}
        
        /* Navigation */
        .nav-tabs {{
            display: flex;
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            border-radius: 10px;
            padding: 10px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            gap: 5px;
        }}
        
        .nav-tab {{
            padding: 12px 20px;
            background: transparent;
            color: var(--light);
            border: none;
            border-radius: 8px;
            cursor: pointer;
            transition: all 0.3s ease;
            font-weight: 500;
        }}
        
        .nav-tab:hover {{
            background: rgba(255, 255, 255, 0.2);
        }}
        
        .nav-tab.active {{
            background: var(--accent);
            color: var(--dark);
        }}
        
        /* Content Sections */
        .content-section {{
            display: none;
            background: rgba(255, 255, 255, 0.95);
            color: var(--dark);
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.2);
            margin-bottom: 20px;
        }}
        
        .content-section.active {{
            display: block;
        }}
        
        .section-title {{
            color: var(--primary);
            border-bottom: 3px solid var(--secondary);
            padding-bottom: 15px;
            margin-bottom: 25px;
            font-size: 28px;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        /* Statistics Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, var(--secondary), var(--primary));
            color: white;
            padding: 25px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        
        .stat-number {{
            font-size: 32px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        
        /* Lineage Grid */
        .lineage-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
            gap: 25px;
        }}
        
        .lineage-card {{
            background: white;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
            overflow: hidden;
            transition: transform 0.3s ease, box-shadow 0.3s ease;
            border-left: 5px solid var(--secondary);
        }}
        
        .lineage-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.2);
        }}
        
        .lineage-header {{
            background: linear-gradient(135deg, var(--primary), var(--secondary));
            color: white;
            padding: 20px;
        }}
        
        .lineage-name {{
            font-size: 24px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .lineage-subtitle {{
            font-size: 16px;
            opacity: 0.9;
            margin-bottom: 10px;
        }}
        
        .lineage-badges {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        
        .badge {{
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: bold;
        }}
        
        .badge-healthcare {{ background: var(--danger); }}
        .badge-community {{ background: var(--success); }}
        .badge-livestock {{ background: var(--warning); color: var(--dark); }}
        .badge-high {{ background: var(--danger); }}
        .badge-medium {{ background: var(--warning); color: var(--dark); }}
        .badge-low {{ background: var(--success); }}
        
        .lineage-content {{
            padding: 20px;
        }}
        
        .info-section {{
            margin-bottom: 15px;
        }}
        
        .info-title {{
            font-weight: bold;
            color: var(--primary);
            margin-bottom: 8px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        
        .gene-list {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-top: 5px;
        }}
        
        .gene-tag {{
            background: var(--secondary);
            color: white;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-family: 'Courier New', monospace;
        }}
        
        .references {{
            font-size: 12px;
            color: #666;
            margin-top: 15px;
            border-top: 1px solid #eee;
            padding-top: 10px;
        }}
        
        /* Specialized Profiles */
        .profile-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: 20px;
        }}
        
        .profile-card {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
            border-left: 4px solid var(--accent);
        }}
        
        /* Geographic Table */
        .geo-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
        }}
        
        .geo-table th {{
            background: var(--primary);
            color: white;
            padding: 15px;
            text-align: left;
        }}
        
        .geo-table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        
        .geo-table tr:hover {{
            background: #f8fafc;
        }}
        
        /* Search and Filter */
        .search-box {{
            background: rgba(255, 255, 255, 0.95);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
        }}
        
        .search-input {{
            width: 100%;
            padding: 15px;
            border: 2px solid var(--secondary);
            border-radius: 8px;
            font-size: 16px;
            margin-bottom: 15px;
        }}
        
        .filter-buttons {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
        }}
        
        .filter-btn {{
            padding: 10px 20px;
            background: var(--secondary);
            color: white;
            border: none;
            border-radius: 20px;
            cursor: pointer;
            transition: background 0.3s ease;
        }}
        
        .filter-btn:hover {{
            background: var(--primary);
        }}
        
        .filter-btn.active {{
            background: var(--accent);
            color: var(--dark);
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            margin-top: 40px;
            padding: 30px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 15px;
        }}
        
        .authorship {{
            margin-top: 20px;
            padding: 20px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 10px;
        }}
        
        @media (max-width: 768px) {{
            .lineage-grid {{
                grid-template-columns: 1fr;
            }}
            .stats-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
            .nav-tabs {{
                flex-direction: column;
            }}
            .ascii-art {{
                font-size: 10px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
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
        </div>

        <!-- Navigation Tabs -->
        <div class="nav-tabs">
            <button class="nav-tab active" data-tab="overview">
                <i class="fas fa-chart-bar"></i> Overview
            </button>
            <button class="nav-tab" data-tab="lineages">
                <i class="fas fa-dna"></i> All Lineages ({total_lineages})
            </button>
            <button class="nav-tab" data-tab="specialized">
                <i class="fas fa-star"></i> Specialized Profiles
            </button>
            <button class="nav-tab" data-tab="geographic">
                <i class="fas fa-globe-americas"></i> Geographic Distribution
            </button>
            <button class="nav-tab" data-tab="about">
                <i class="fas fa-info-circle"></i> About
            </button>
        </div>

        <!-- Search Box -->
        <div class="search-box" id="searchSection">
            <input type="text" class="search-input" id="searchInput" 
                   placeholder="🔍 Search lineages by CC, name, type, genes, or characteristics...">
            <div class="filter-buttons">
                <button class="filter-btn active" data-filter="all">All</button>
                <button class="filter-btn" data-filter="healthcare">HA-MRSA</button>
                <button class="filter-btn" data-filter="community">CA-MRSA</button>
                <button class="filter-btn" data-filter="livestock">LA-MRSA</button>
                <button class="filter-btn" data-filter="high">High Risk</button>
            </div>
        </div>

        <!-- Overview Section -->
        <div class="content-section active" id="overview">
            <h2 class="section-title">
                <i class="fas fa-chart-bar"></i> Database Overview
            </h2>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-number">{total_lineages}</div>
                    <div class="stat-label">Total Lineages</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['healthcare']}</div>
                    <div class="stat-label">HA-MRSA Lineages</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['community']}</div>
                    <div class="stat-label">CA-MRSA Lineages</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['livestock']}</div>
                    <div class="stat-label">LA-MRSA Lineages</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['high_risk']}</div>
                    <div class="stat-label">High Risk Lineages</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{stats['pvl_positive']}</div>
                    <div class="stat-label">PVL-Positive</div>
                </div>
            </div>

            <div style="background: #fff3cd; color: #856404; padding: 20px; border-radius: 8px; margin: 20px 0;">
                <h3 style="color: #856404; margin-bottom: 10px;">
                    <i class="fas fa-exclamation-triangle"></i> Important Usage Notes
                </h3>
                <ul style="margin-left: 20px;">
                    <li>This reference database should be used for comparison and education only</li>
                    <li>Always verify lineage assignments with additional methods and recent literature</li>
                    <li>Lineage distributions and characteristics can vary geographically and change over time</li>
                    <li>Gene presence varies within lineages - not all strains carry all listed virulence factors</li>
                    <li>Refer to PubMLST for current distributions and additional sequence types</li>
                    <li>Note that not all references have been attached and verify any classifications from literature</li>
                </ul>
            </div>
        </div>

        <!-- Lineages Section -->
        <div class="content-section" id="lineages">
            <h2 class="section-title">
                <i class="fas fa-dna"></i> Complete Lineage Database
            </h2>
            <div class="lineage-grid" id="lineageGrid">
'''

    # Generate lineage cards
    for cc, info in sorted(lineages.items()):
        # Determine badge classes
        type_class = ""
        type_icon = ""
        if "Healthcare" in info.get("type", ""):
            type_class = "badge-healthcare"
            type_icon = "🏥"
        elif "Community" in info.get("type", ""):
            type_class = "badge-community" 
            type_icon = "🏘️"
        elif "Livestock" in info.get("type", ""):
            type_class = "badge-livestock"
            type_icon = "🐄"
        
        risk_class = f"badge-{info.get('risk_level', 'medium').lower()}"
        
        html_content += f'''
                <div class="lineage-card" data-type="{info.get('type', '').lower()}" data-risk="{info.get('risk_level', 'MEDIUM').lower()}">
                    <div class="lineage-header">
                        <div class="lineage-name">{cc}</div>
                        <div class="lineage-subtitle">{info.get('primary_name', 'Unknown')}</div>
                        <div class="lineage-badges">
                            <span class="badge {type_class}">{type_icon} {info.get('type', 'Unknown')}</span>
                            <span class="badge {risk_class}">Risk: {info.get('risk_level', 'Unknown')}</span>
                        </div>
                    </div>
                    
                    <div class="lineage-content">
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-list-ol"></i> Sequence Types</div>
                            {', '.join(map(str, info.get('sequence_types', [])))}
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-tag"></i> Common spa Types</div>
                            {', '.join(info.get('common_spa_types', []))}
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-cog"></i> SCCmec Types</div>
                            {', '.join(info.get('sccmec_types', []))}
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-map-marker-alt"></i> Geographic Distribution</div>
                            {info.get('geographic_distribution', {}).get('notes', '')}
                            <div style="font-size: 12px; margin-top: 5px;">
                                <strong>Regions:</strong> {', '.join(info.get('geographic_distribution', {}).get('regions', []))}
                            </div>
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-biohazard"></i> Virulence Profile</div>
                            <div class="gene-list">
        '''
        
        # Add all virulence factors
        vp = info.get('virulence_profile', {})
        for category, genes in vp.items():
            if isinstance(genes, list):
                for gene in genes:
                    html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                            </div>
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-shield-alt"></i> Resistance Profile</div>
                            <div><strong>Classes:</strong> {', '.join(info.get('resistance_profile', {}).get('antimicrobial_classes', []))}</div>
                            <div class="gene-list" style="margin-top: 8px;">
        '''
        
        # Add resistance genes
        for gene in info.get('resistance_profile', {}).get('common_genes', []):
            html_content += f'<span class="gene-tag">{gene}</span>'
        
        html_content += f'''
                            </div>
                        </div>
                        
                        <div class="info-section">
                            <div class="info-title"><i class="fas fa-stethoscope"></i> Clinical Significance</div>
                            {info.get('clinical_significance', '')}
                        </div>
                        
                        <div class="references">
                            <strong>Key References:</strong> {', '.join(info.get('key_references', []))}
                        </div>
                    </div>
                </div>
        '''
    
    html_content += '''            </div>
        </div>

        <!-- Specialized Profiles Section -->
        <div class="content-section" id="specialized">
            <h2 class="section-title">
                <i class="fas fa-star"></i> Specialized Lineage Profiles
            </h2>
'''
    
    # Generate specialized profiles
    for profile_type, profiles in specialized.items():
        html_content += f'''
            <h3 style="color: var(--primary); margin: 25px 0 15px 0; padding-bottom: 10px; border-bottom: 2px solid var(--secondary);">
                {profile_type.replace('_', ' ').title()}
            </h3>
            <div class="profile-grid">
        '''
        
        for profile_name, profile_data in profiles.items():
            html_content += f'''
                <div class="profile-card">
                    <h4 style="color: var(--primary); margin-bottom: 10px;">{profile_name}</h4>
                    <div style="margin-bottom: 8px;"><strong>ST:</strong> {profile_data.get('st', 'N/A')}</div>
                    <div style="margin-bottom: 8px;"><strong>SCCmec:</strong> {profile_data.get('sccmec', 'N/A')}</div>
                    <div style="margin-bottom: 8px;"><strong>spa:</strong> {profile_data.get('spa', 'N/A')}</div>
                    <div style="margin-bottom: 8px;"><strong>Clinical:</strong> {profile_data.get('clinical', 'N/A')}</div>
                    <div style="margin-bottom: 8px;"><strong>Risk:</strong> {profile_data.get('risk', 'N/A')}</div>
                    <div class="gene-list" style="margin-top: 10px;">
            '''
            
            for virulence in profile_data.get('virulence', []):
                html_content += f'<span class="gene-tag">{virulence}</span>'
            
            if profile_data.get('resistance'):
                html_content += '<div style="margin-top: 8px;"><strong>Resistance:</strong> '
                for resistance in profile_data.get('resistance', []):
                    html_content += f'<span class="gene-tag">{resistance}</span>'
                html_content += '</div>'
            
            html_content += '''
                    </div>
                </div>
            '''
        
        html_content += '''            </div>
        '''
    
    html_content += '''        </div>

        <!-- Geographic Distribution Section -->
        <div class="content-section" id="geographic">
            <h2 class="section-title">
                <i class="fas fa-globe-americas"></i> Global Geographic Distribution
            </h2>
            
            <table class="geo-table">
                <thead>
                    <tr>
                        <th>Region</th>
                        <th>Lineage/Clone</th>
                        <th>Clonal Complex</th>
                        <th>Prevalence</th>
                        <th>Setting</th>
                    </tr>
                </thead>
                <tbody>
'''
    
    # Generate geographic table
    for region, lineages_data in geographic.items():
        for clone_name, clone_data in lineages_data.items():
            html_content += f'''
                    <tr>
                        <td><strong>{region.replace('_', ' ').title()}</strong></td>
                        <td>{clone_name}</td>
                        <td>{clone_data.get('cc', 'N/A')}</td>
                        <td>{clone_data.get('prevalence', 'N/A')}</td>
                        <td>{clone_data.get('setting', 'N/A')}</td>
                    </tr>
            '''
    
    html_content += '''                </tbody>
            </table>
        </div>

        <!-- About Section -->
        <div class="content-section" id="about">
            <h2 class="section-title">
                <i class="fas fa-info-circle"></i> About StaphScope Lineage Database
            </h2>
            
            <div style="background: white; padding: 25px; border-radius: 10px; margin-bottom: 20px;">
                <h3 style="color: var(--primary); margin-bottom: 15px;">Database Information</h3>
                <p><strong>Version:</strong> 2025</p>
                <p><strong>Last Updated:</strong> 20/10/2025</p>
                <p><strong>Total Lineages:</strong> 26</p>
                <p><strong>Coverage:</strong> Global S. aureus lineages including HA-MRSA, CA-MRSA, and LA-MRSA</p>
                
                <h4 style="color: var(--primary); margin: 20px 0 10px 0;">Database Structure</h4>
                <ul style="margin-left: 20px;">
                    <li><strong>Lineage Database:</strong> clonal complexes with detailed metadata</li>
                    <li><strong>Specialized Profiles:</strong> categories of specialized lineages</li>
                    <li><strong>Geographic Data:</strong> regions with lineage distributions</li>
                </ul>
            </div>
            
            <div class="authorship">
                <h3 style="color: var(--accent); margin-bottom: 15px;">Technical Support & Inquiries</h3>
                <p><strong>Author:</strong> Brown Beckley</p>
                <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" style="color: var(--accent);">bbeckley-hub</a></p>
                <p><strong>Email:</strong> <a href="mailto:brownbeckley94@gmail.com" style="color: var(--accent);">brownbeckley94@gmail.com</a></p>
                <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
                
                <div style="margin-top: 20px; padding: 15px; background: rgba(0,0,0,0.1); border-radius: 8px;">
                    <p style="font-size: 14px;"><strong>Note:</strong> This database is for reference and educational purposes. Always verify lineage assignments with current literature, Pubmlst and additional molecular methods.</p>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Tab navigation
        const navTabs = document.querySelectorAll('.nav-tab');
        const contentSections = document.querySelectorAll('.content-section');
        
        navTabs.forEach(tab => {{
            tab.addEventListener('click', () => {{
                const targetTab = tab.dataset.tab;
                
                // Update active tab
                navTabs.forEach(t => t.classList.remove('active'));
                tab.classList.add('active');
                
                // Show target section
                contentSections.forEach(section => {{
                    section.classList.remove('active');
                    if (section.id === targetTab) {{
                        section.classList.add('active');
                    }}
                }});
                
                // Show/hide search based on section
                const searchSection = document.getElementById('searchSection');
                if (targetTab === 'lineages') {{
                    searchSection.style.display = 'block';
                }} else {{
                    searchSection.style.display = 'none';
                }}
            }});
        }});
        
        // Search and filter functionality
        const searchInput = document.getElementById('searchInput');
        const filterButtons = document.querySelectorAll('.filter-btn');
        const lineageCards = document.querySelectorAll('.lineage-card');
        
        function filterLineages() {{
            const searchTerm = searchInput.value.toLowerCase();
            const activeFilter = document.querySelector('.filter-btn.active').dataset.filter;
            
            lineageCards.forEach(card => {{
                const cardText = card.textContent.toLowerCase();
                const cardType = card.dataset.type;
                const cardRisk = card.dataset.risk;
                
                const matchesSearch = searchTerm === '' || cardText.includes(searchTerm);
                const matchesFilter = 
                    activeFilter === 'all' ||
                    (activeFilter === 'healthcare' && cardType.includes('healthcare')) ||
                    (activeFilter === 'community' && cardType.includes('community')) ||
                    (activeFilter === 'livestock' && cardType.includes('livestock')) ||
                    (activeFilter === 'high' && cardRisk === 'high');
                
                card.style.display = (matchesSearch && matchesFilter) ? 'block' : 'none';
            }});
        }}
        
        searchInput.addEventListener('input', filterLineages);
        
        filterButtons.forEach(btn => {{
            btn.addEventListener('click', () => {{
                filterButtons.forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                filterLineages();
            }});
        }});
        
        // Initialize
        filterLineages();
    </script>
</body>
</html>
'''
    
    # Save HTML file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Comprehensive HTML reference generated: {output_path}")
    return output_path

# Quick test function
def test_html_generation():
    """Test the HTML generation"""
    output_path = generate_comprehensive_lineage_html()
    print(f"📁 File saved at: {os.path.abspath(output_path)}")
    print(f"📊 Open the file in your browser to view the reference database")
    return output_path

if __name__ == "__main__":
    test_html_generation()
