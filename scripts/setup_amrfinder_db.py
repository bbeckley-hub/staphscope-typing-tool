#!/usr/bin/env python3
"""
Setup AMRFinderPlus database automatically
"""
import subprocess
import sys
import os

def setup_amrfinder_database():
    """Setup AMRFinderPlus database"""
    print("🔧 Setting up AMRFinderPlus database...")
    
    try:
        # Update AMRFinderPlus database
        result = subprocess.run(
            ["amrfinder", "--update"],
            capture_output=True,
            text=True,
            check=True
        )
        print("✅ AMRFinderPlus database update completed")
        
        # Verify the database
        result = subprocess.run(
            ["amrfinder", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"✅ AMRFinderPlus version: {result.stdout.strip()}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error setting up AMRFinderPlus database: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("❌ Error: amrfinder not found in PATH")
        return False

if __name__ == "__main__":
    success = setup_amrfinder_database()
    sys.exit(0 if success else 1)
