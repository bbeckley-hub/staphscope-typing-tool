#!/usr/bin/env python3
"""
StaphScope Banner Module
Beautiful ASCII art and scientific quotes for terminal display
Author: Brown Beckley <brownbeckley94@gmail.com>
"""

import random
from datetime import datetime
import sys
import time
import textwrap

class StaphScopeBanner:
    """StaphScope Banner Display with Scientific Quotes"""
    
    def __init__(self):
        self.banner_art = self._get_banner_art()
        self.quotes = self._get_scientific_quotes()
        self.version = "v1.0."
        self.author_info = self._get_author_info()
    
    def _get_banner_art(self):
        """Return the main StaphScope ASCII art"""
        return r"""
    ███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
    ██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
    ███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
    ╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
    ███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
    ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
    
        🧫 Advanced Staphylococcus aureus Typing & Lineage Analysis Platform
    """
    
    def _get_scientific_quotes(self):
        """Return collection of scientific quotes about microbiology and genomics"""
        return [
            {
                "quote": "The important thing in science is not so much to obtain new facts as to discover new ways of thinking about them.",
                "author": "William Lawrence Bragg"
            },
            {
                "quote": "In the fields of observation chance favors only the prepared mind.",
                "author": "Louis Pasteur"
            },
            {
                "quote": "Nothing in life is to be feared, it is only to be understood. Now is the time to understand more, so that we may fear less.",
                "author": "Marie Curie"
            },
            {
                "quote": "The genome is the book of life, and we are learning to read it with increasing clarity.",
                "author": "Francis Collins"
            },
            {
                "quote": "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.",
                "author": "Stephen Hawking"
            },
            {
                "quote": "We are all star-stuff contemplating the stars.",
                "author": "Carl Sagan"
            },
            {
                "quote": "The microbe is nothing; the terrain is everything.",
                "author": "Louis Pasteur"
            },
            {
                "quote": "DNA is like a computer program but far, far more advanced than any software ever created.",
                "author": "Bill Gates"
            },
            {
                "quote": "The art of medicine consists of amusing the patient while nature cures the disease.",
                "author": "Voltaire"
            },
            {
                "quote": "We are just an advanced breed of monkeys on a minor planet of a very average star. But we can understand the Universe. That makes us something very special.",
                "author": "Stephen Hawking"
            }
        ]
    
    def _get_author_info(self):
        """Return author information"""
        return {
            "name": "Brown Beckley",
            "github": "bbeckley-hub",
            "email": "brownbeckley94@gmail.com",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "license": "MIT"
        }
    
    def _get_colors(self):
        """Define color codes for terminal output"""
        class Colors:
            RED = '\033[91m'
            GREEN = '\033[92m'
            YELLOW = '\033[93m'
            BLUE = '\033[94m'
            MAGENTA = '\033[95m'
            CYAN = '\033[96m'
            WHITE = '\033[97m'
            BOLD = '\033[1m'
            UNDERLINE = '\033[4m'
            END = '\033[0m'
        return Colors
    
    def display_banner(self, show_quote=True, show_author=True):
        """Display the main StaphScope banner"""
        C = self._get_colors()
        
        print(f"{C.CYAN}{C.BOLD}{self.banner_art}{C.END}")
        print(f"{C.YELLOW}{C.BOLD}         Version: {self.version}{C.END}")
        
        if show_quote:
            quote = random.choice(self.quotes)
            # Format quote without box
            print(f"{C.GREEN}╔══════════════════════════════════════════════════════════════════════════════╗{C.END}")
            print(f"{C.GREEN}║ {C.CYAN}💡 SCIENTIFIC QUOTE:{C.END}")
            
            # Wrap long quotes to fit within the box
            wrapped_quote = textwrap.fill(quote['quote'], width=72)
            quote_lines = wrapped_quote.split('\n')
            
            for line in quote_lines:
                print(f"{C.GREEN}║ {C.WHITE}{line:<72}{C.GREEN}║{C.END}")
            
            print(f"{C.GREEN}║ {C.YELLOW}   — {quote['author']:<65} {C.GREEN}║{C.END}")
            print(f"{C.GREEN}╚══════════════════════════════════════════════════════════════════════════════╝{C.END}")
            print()
        
        if show_author:
            print(f"{C.GREEN}╔══════════════════════════════════════════════════════════════════════════════╗{C.END}")
            print(f"{C.GREEN}║ {C.MAGENTA}👨‍💻 AUTHOR INFORMATION:{C.END}")
            print(f"{C.GREEN}║ {C.MAGENTA}   Name: {self.author_info['name']:<60} {C.GREEN}║{C.END}")
            print(f"{C.GREEN}║ {C.BLUE}   GitHub: {self.author_info['github']:<58} {C.GREEN}║{C.END}")
            print(f"{C.GREEN}║ {C.CYAN}   Email: {self.author_info['email']:<58} {C.GREEN}║{C.END}")
            
            # Handle long affiliation by splitting if needed
            affiliation = self.author_info['affiliation']
            if len(affiliation) > 60:
                # Split affiliation into two lines if too long
                affil_parts = textwrap.wrap(affiliation, width=60)
                print(f"{C.GREEN}║ {C.WHITE}   Affiliation: {affil_parts[0]:<53} {C.GREEN}║{C.END}")
                for part in affil_parts[1:]:
                    print(f"{C.GREEN}║ {C.WHITE}               {part:<53} {C.GREEN}║{C.END}")
            else:
                print(f"{C.GREEN}║ {C.WHITE}   Affiliation: {affiliation:<53} {C.GREEN}║{C.END}")
            
            print(f"{C.GREEN}║ {C.YELLOW}   License: {self.author_info['license']:<57} {C.GREEN}║{C.END}")
            print(f"{C.GREEN}╚══════════════════════════════════════════════════════════════════════════════╝{C.END}")
        print()
    
    def display_startup_sequence(self):
        """Display animated startup sequence"""
        C = self._get_colors()
        
        print(f"{C.CYAN}{C.BOLD}🚀 Initializing StaphScope Analysis Platform...{C.END}")
        time.sleep(0.5)
        
        steps = [
            "Loading genomic databases...",
            "Initializing MLST analysis engine...",
            "Configuring lineage prediction algorithms...",
            "Setting up resistance gene detection...",
            "Preparing virulence factor analysis...",
            "Setting up spatyping prediction...",
            "Enabling Sccemecfinder analysis...",
            "Optimizing multi-threading capabilities...",
        ]
        
        for i, step in enumerate(steps, 1):
            print(f"{C.YELLOW}[{i}/{len(steps)}] {step}{C.END}")
            time.sleep(0.3)
        
        print(f"{C.GREEN}{C.BOLD}✅ StaphScope ready for analysis!{C.END}")
        print()
    
    def display_footer(self, analysis_time=None, samples_processed=0):
        """Display analysis completion footer"""
        C = self._get_colors()
        
        print()
        print(f"{C.MAGENTA}{C.BOLD}╔══════════════════════════════════════════════════════════════════════════════╗{C.END}")
        print(f"{C.MAGENTA}{C.BOLD}║                        ANALYSIS COMPLETE - STAPHSCOPE                        ║{C.END}")
        print(f"{C.MAGENTA}{C.BOLD}╠══════════════════════════════════════════════════════════════════════════════╣{C.END}")
        
        if analysis_time:
            print(f"{C.MAGENTA}║ {C.CYAN}⏱️  Analysis Duration: {analysis_time:<50} {C.MAGENTA}║{C.END}")
        
        if samples_processed > 0:
            print(f"{C.MAGENTA}║ {C.GREEN}🧫 Samples Processed: {samples_processed:<49} {C.MAGENTA}║{C.END}")
        
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{C.MAGENTA}║ {C.YELLOW}📅 Completion Time: {current_time:<48} {C.MAGENTA}║{C.END}")
        print(f"{C.MAGENTA}╠══════════════════════════════════════════════════════════════════════════════╣{C.END}")
        print(f"{C.MAGENTA}║ {C.WHITE}TECHNICAL SUPPORT & INQUIRIES:                                        {C.MAGENTA}║{C.END}")
        print(f"{C.MAGENTA}║ {C.CYAN}   Author: {self.author_info['name']:<55} {C.MAGENTA}║{C.END}")
        print(f"{C.MAGENTA}║ {C.BLUE}   GitHub: {self.author_info['github']:<56} {C.MAGENTA}║{C.END}")
        print(f"{C.MAGENTA}║ {C.GREEN}   Email: {self.author_info['email']:<56} {C.MAGENTA}║{C.END}")
        
        # Handle long affiliation in footer
        affiliation = self.author_info['affiliation']
        if len(affiliation) > 60:
            affil_parts = textwrap.wrap(affiliation, width=60)
            print(f"{C.MAGENTA}║ {C.WHITE}   Affiliation: {affil_parts[0]:<52} {C.MAGENTA}║{C.END}")
            for part in affil_parts[1:]:
                print(f"{C.MAGENTA}║ {C.WHITE}               {part:<52} {C.MAGENTA}║{C.END}")
        else:
            print(f"{C.MAGENTA}║ {C.WHITE}   Affiliation: {affiliation:<52} {C.MAGENTA}║{C.END}")
        
        print(f"{C.MAGENTA}╚══════════════════════════════════════════════════════════════════════════════╝{C.END}")
        print()
    
    def display_module_header(self, module_name, description=""):
        """Display header for specific analysis modules"""
        C = self._get_colors()
        
        print()
        print(f"{C.BLUE}{C.BOLD}╔══════════════════════════════════════════════════════════════════════════════╗{C.END}")
        print(f"{C.BLUE}{C.BOLD}║ 🧬 {module_name.upper():<70} ║{C.END}")
        if description:
            print(f"{C.BLUE}{C.BOLD}║   {description:<68} ║{C.END}")
        print(f"{C.BLUE}{C.BOLD}╚══════════════════════════════════════════════════════════════════════════════╝{C.END}")
        print()
    
    def display_warning(self, message):
        """Display warning message"""
        C = self._get_colors()
        print(f"{C.YELLOW}⚠️  WARNING: {message}{C.END}")
    
    def display_error(self, message):
        """Display error message"""
        C = self._get_colors()
        print(f"{C.RED}❌ ERROR: {message}{C.END}")
    
    def display_success(self, message):
        """Display success message"""
        C = self._get_colors()
        print(f"{C.GREEN}✅ SUCCESS: {message}{C.END}")
    
    def display_info(self, message):
        """Display info message"""
        C = self._get_colors()
        print(f"{C.CYAN}💡 INFO: {message}{C.END}")
    
    def display_progress_bar(self, iteration, total, prefix='', suffix='', length=50, fill='█'):
        """Display progress bar"""
        C = self._get_colors()
        percent = ("{0:.1f}").format(100 * (iteration / float(total)))
        filled_length = int(length * iteration // total)
        bar = fill * filled_length + '-' * (length - filled_length)
        print(f'\r{C.CYAN}{prefix} |{bar}| {percent}% {suffix}{C.END}', end='\r')
        if iteration == total:
            print()

def main():
    """Test the banner display"""
    banner = StaphScopeBanner()
    
    # Display full banner with startup sequence
    banner.display_startup_sequence()
    banner.display_banner(show_quote=True, show_author=True)
    
    # Test module headers
    banner.display_module_header("MLST Analysis", "Multi-Locus Sequence Typing for S. aureus")
    banner.display_success("MLST analysis completed successfully!")
    
    banner.display_module_header("Lineage Prediction", "Clonal complex and epidemic clone identification")
    banner.display_info("Processing 5 samples with multi-threading")
    
    # Test progress bar
    print("Progress demonstration:")
    for i in range(101):
        banner.display_progress_bar(i, 100, prefix='Progress:', suffix='Complete', length=40)
        time.sleep(0.02)
    print()
    
    # Display footer
    banner.display_footer(analysis_time="2 minutes, 15 seconds", samples_processed=5)

if __name__ == "__main__":
    main()
