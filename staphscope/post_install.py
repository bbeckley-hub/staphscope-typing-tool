#!/usr/bin/env python3
"""
Post-installation script for StaphScope
Runs necessary database setup commands
"""
import os
import subprocess
import sys
import time

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n🔧 {description}...")
    print(f"Running: {cmd}")
    
    try:
        # Run the command
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"✅ {description} completed successfully")
            if result.stdout:
                print(f"Output: {result.stdout}")
        else:
            print(f"⚠️ {description} had issues (return code: {result.returncode})")
            if result.stderr:
                print(f"Error: {result.stderr}")
            # Don't fail the installation for database setup issues
            return False
    except Exception as e:
        print(f"❌ Error running {description}: {str(e)}")
        return False
    
    return True

def main():
    print("🚀 StaphScope Post-Installation Setup")
    print("=" * 50)
    
    # Get current directory (where package is installed)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"Package installed in: {current_dir}")
    
    # 1. Setup AMRFinder database
    run_command("amrfinder -u", "Updating AMRFinder database")
    
    # 2. Setup Abricate databases
    run_command("abricate --setupdb", "Setting up Abricate databases")
    
    # 3. Add current directory to PYTHONPATH
    python_path_addition = f'export PYTHONPATH=$PYTHONPATH:"{current_dir}"'
    
    # Create a setup script for users
    setup_script = f"""#!/bin/bash
# StaphScope environment setup
echo "Setting up StaphScope environment..."
export PYTHONPATH=$PYTHONPATH:"{current_dir}"
echo "Added {current_dir} to PYTHONPATH"
echo "StaphScope setup complete!"
"""
    
    setup_script_path = os.path.join(current_dir, "setup_staphscope_env.sh")
    with open(setup_script_path, 'w') as f:
        f.write(setup_script)
    
    os.chmod(setup_script_path, 0o755)
    print(f"✅ Created environment setup script: {setup_script_path}")
    
    # Print instructions for users
    print("\n📋 Post-Installation Instructions:")
    print("1. Run the following command to set up your environment:")
    print(f"   source {setup_script_path}")
    print("2. Or manually add to your ~/.bashrc or ~/.zshrc:")
    print(f'   export PYTHONPATH=$PYTHONPATH:"{current_dir}"')
    
    print("\n🎉 StaphScope installation completed!")
    print("You can now run: staphscope --help")

if __name__ == "__main__":
    main()
