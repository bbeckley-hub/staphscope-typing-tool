#!/usr/bin/env python3
"""
Post-installation setup for StaphScope
"""
import subprocess
import sys
import os

def run_post_install():
    """Run all post-installation setup tasks"""
    print("🚀 Running StaphScope post-installation setup...")
    print("This may take several minutes as we download and setup all databases...")
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    results = []
    
    # Run abricate database setup (ALL databases)
    abricate_script = os.path.join(script_dir, "setup_abricate_db.py")
    print("\n" + "="*60)
    print("Setting up ALL Abricate databases...")
    result1 = subprocess.run([sys.executable, abricate_script])
    results.append(result1.returncode)
    
    # Run AMRFinderPlus database setup
    amrfinder_script = os.path.join(script_dir, "setup_amrfinder_db.py")
    print("\n" + "="*60)
    print("Setting up AMRFinderPlus database...")
    result2 = subprocess.run([sys.executable, amrfinder_script])
    results.append(result2.returncode)
    
    print("\n" + "="*60)
    success_count = sum(1 for r in results if r == 0)
    
    if success_count == len(results):
        print("🎉 StaphScope installation completed SUCCESSFULLY!")
        print("💡 All databases are ready. You can now run: staphscope -i input.fasta -o results")
    elif success_count >= 1:
        print("✅ StaphScope installation completed with partial success.")
        print(f"💡 {success_count}/{len(results)} database setups completed successfully.")
        print("💡 You can now run: staphscope -i input.fasta -o results")
    else:
        print("⚠️  StaphScope installation completed with issues.")
        print("💡 You may need to manually setup some databases.")
    
    return success_count > 0  # Success if at least one DB setup worked

if __name__ == "__main__":
    success = run_post_install()
    sys.exit(0 if success else 1)
