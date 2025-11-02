#!/usr/bin/env python3
"""
Setup ALL abricate databases automatically
"""
import subprocess
import sys
import os

def setup_abricate_databases():
    """Setup ALL abricate databases"""
    print("🔧 Setting up ALL abricate databases...")
    
    # ALL abricate databases
    databases = [
        "resfinder", "card", "ncbi", "argannot", "megares", 
        "ecoh", "plasmidfinder", "ecoli_vf", "vfdb"
    ]
    
    successful_dbs = []
    failed_dbs = []
    
    for db in databases:
        print(f"📦 Setting up {db} database...")
        try:
            result = subprocess.run(
                ["abricate", "--setupdb", db],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"✅ {db} database setup completed")
            successful_dbs.append(db)
        except subprocess.CalledProcessError as e:
            print(f"⚠️  Could not setup {db} database: {e}")
            failed_dbs.append(db)
        except FileNotFoundError:
            print("❌ Error: abricate not found in PATH")
            return False
    
    # Verify setup
    try:
        result = subprocess.run(["abricate", "--list"], capture_output=True, text=True, check=True)
        print("\n✅ Abricate databases setup summary:")
        print(f"   Successfully setup: {len(successful_dbs)} databases")
        print(f"   Failed to setup: {len(failed_dbs)} databases")
        if successful_dbs:
            print(f"   Available databases: {', '.join(successful_dbs)}")
        if failed_dbs:
            print(f"   Failed databases: {', '.join(failed_dbs)}")
        
        print("\nFull database list:")
        print(result.stdout)
        return len(successful_dbs) > 0  # Success if at least one DB works
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error verifying abricate setup: {e}")
        return False

if __name__ == "__main__":
    success = setup_abricate_databases()
    sys.exit(0 if success else 1)
