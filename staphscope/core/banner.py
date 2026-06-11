#!/usr/bin/env python3
"""
StaphScope Banner Module - FIXED VERSION WITH 8-STEP STARTUP
Beautiful ASCII art and scientific quotes for terminal display
Author: Brown Beckley <brownbeckley94@gmail.com>
Date: 2025
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Biochemistry
version 1.2.3(2026-06-10)
"""

import random
from datetime import datetime
import sys
import time
import textwrap
import os
import shutil

class StaphScopeBanner:
    """StaphScope Banner Display with Scientific Quotes"""
    
    def __init__(self):
        self.banner_art = self._get_banner_art()
        self.quotes = self._get_scientific_quotes()
        self.version = "v1.2.3"
        self.author_info = self._get_author_info()
        self.terminal_width = self._get_terminal_width()
    
    def _get_terminal_width(self):
        """Get terminal width, default to 88 if cannot determine"""
        try:
            return min(100, shutil.get_terminal_size().columns - 2)
        except:
            return 88
    
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
            DIM = '\033[2m'
            END = '\033[0m'
            # Extended colors
            LIGHT_BLUE = '\033[38;5;117m'
            LIGHT_GREEN = '\033[38;5;120m'
            LIGHT_YELLOW = '\033[38;5;228m'
            ORANGE = '\033[38;5;214m'
            PURPLE = '\033[38;5;141m'
            PINK = '\033[38;5;213m'
        return Colors
    
    def display_banner(self, show_quote=True, show_author=True):
        """Display the main StaphScope banner"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Main banner
        print(f"\n{C.CYAN}{C.BOLD}{self.banner_art}{C.END}")
        
        # Top decoration
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        
        # Version and date
        version_text = f"🔬 Version: {self.version}  |  {datetime.now().strftime('%Y-%m-%d')}  |  StaphScope Platform 🧬"
        print(f"{C.YELLOW}{C.BOLD}{version_text.center(width)}{C.END}")
        
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        if show_quote:
            quote = random.choice(self.quotes)
            
            # Quote section header
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{C.BOLD}{'💡 SCIENTIFIC INSPIRATION'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            # Quote text with wrapping
            wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
            quote_lines = wrapper.wrap(f'"{quote["quote"]}"')
            
            for line in quote_lines:
                print(f"{C.WHITE}{line}{C.END}")
            
            print()
            author_text = f"— {quote['author']}"
            print(f"{C.YELLOW}{C.BOLD}{author_text.rjust(width-2)}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
        
        if show_author:
            # Author section header
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{C.BOLD}{'👨‍💻 DEVELOPER & CONTACT INFORMATION'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            # Author details in columns
            print(f"{C.LIGHT_BLUE}  👤 Name:{C.END}        {C.WHITE}{self.author_info['name']}{C.END}")
            print(f"{C.LIGHT_BLUE}  🐙 GitHub:{C.END}      {C.WHITE}{self.author_info['github']}{C.END}")
            print(f"{C.LIGHT_BLUE}  📧 Email:{C.END}       {C.WHITE}{self.author_info['email']}{C.END}")
            print(f"{C.LIGHT_BLUE}  🏛️  Affiliation:{C.END} {C.WHITE}{self.author_info['affiliation']}{C.END}")
            print(f"{C.LIGHT_BLUE}  📜 License:{C.END}     {C.WHITE}{self.author_info['license']}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
    
    def display_startup_sequence(self):
        """Display animated startup sequence with progress - USING THE 8-STEP VERSION YOU LIKED"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Header
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.CYAN}{'█' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'🚀 INITIALIZING STAPHSCOPE ANALYSIS PLATFORM 🚀'.center(width)}{C.END}")
        print(f"{C.CYAN}{'█' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}\n")
        time.sleep(0.3)
        
        steps = [
            ("🗄️", "Loading genomic databases", C.BLUE),
            ("🧬", "Initializing MLST analysis engine", C.CYAN),
            ("🌳", "Configuring lineage prediction algorithms", C.GREEN),
            ("🛡️", "Setting up resistance gene detection", C.YELLOW),
            ("⚔️", "Preparing virulence factor analysis", C.ORANGE),
            ("🔬", "Setting up spa-typing prediction", C.MAGENTA),
            ("🔍", "Enabling SCCmec finder analysis", C.LIGHT_BLUE),
            ("⚡", "Optimizing multi-threading capabilities", C.LIGHT_GREEN),
        ]
        
        for i, (icon, step, color) in enumerate(steps, 1):
            print(f"{color}[{i}/{len(steps)}] {icon}  {step}...{C.END}")
            # Progress animation - Using the exact animation from your original code
            for j in range(3):
                print(f"{color}{'  ▓' * (j + 1)}{'░' * (20 - (j + 1) * 3)}{C.END}", end='\r')
                time.sleep(0.08)
            print(f"{color}{'  ▓' * 20} {C.GREEN}✓{C.END}")
            time.sleep(0.1)
        
        # Footer
        print(f"\n{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.GREEN}{'▓' * width}{C.END}")
        print(f"{C.GREEN}{C.BOLD}{'✅ STAPHSCOPE READY FOR ANALYSIS! ✅'.center(width)}{C.END}")
        print(f"{C.GREEN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_GREEN}{'═' * width}{C.END}\n")
        time.sleep(0.3)
    
    def display_footer(self, analysis_time=None, samples_processed=0):
        """Display analysis completion footer"""
        C = self._get_colors()
        width = self.terminal_width
        
        print()
        # Header
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'━' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'🎉 ANALYSIS COMPLETE 🎉'.center(width)}{C.END}")
        print(f"{C.LIGHT_BLUE}{'━' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        # Statistics
        if analysis_time or samples_processed > 0:
            print(f"{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            print(f"{C.YELLOW}{C.BOLD}{'📊 STATISTICS'.center(width)}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            print(f"{C.CYAN}{'▓' * width}{C.END}\n")
            
            if analysis_time:
                print(f"{C.WHITE}  ⏱️  Analysis Duration: {C.GREEN}{C.BOLD}{analysis_time}{C.END}")
            if samples_processed > 0:
                print(f"{C.WHITE}  🧫 Samples Processed:  {C.GREEN}{C.BOLD}{samples_processed}{C.END}")
            
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{C.WHITE}  📅 Completion Time:    {C.CYAN}{C.BOLD}{current_time}{C.END}")
            
            print(f"\n{C.CYAN}{'▓' * width}{C.END}")
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
        
        # Support section
        print(f"{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{C.BOLD}{'📞 SUPPORT & INQUIRIES'.center(width)}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        print(f"{C.WHITE}  👤 Contact:  {C.LIGHT_BLUE}{self.author_info['name']}{C.END}")
        print(f"{C.WHITE}  🐙 GitHub:   {C.LIGHT_BLUE}{self.author_info['github']}{C.END}")
        print(f"{C.WHITE}  📧 Email:    {C.LIGHT_BLUE}{self.author_info['email']}{C.END}")
        
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}\n")
    
    def display_module_header(self, module_name, description=""):
        """Display header for specific analysis modules"""
        C = self._get_colors()
        width = self.terminal_width
        
        print(f"\n{C.CYAN}{'▓' * width}{C.END}")
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        
        # Module title
        title = f"🧬 {module_name.upper()}"
        print(f"{C.YELLOW}{C.BOLD}{title.center(width)}{C.END}")
        
        if description:
            print(f"{C.LIGHT_BLUE}{'─' * width}{C.END}")
            # Wrap long descriptions
            wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
            desc_lines = wrapper.wrap(description)
            for line in desc_lines:
                print(f"{C.WHITE}{line}{C.END}")
        
        print(f"{C.LIGHT_BLUE}{'═' * width}{C.END}")
        print(f"{C.CYAN}{'▓' * width}{C.END}\n")
        
        # Flush output to ensure header displays immediately
        sys.stdout.flush()
    
    def display_section_divider(self, title="", color=None):
        """Display a colorful section divider"""
        C = self._get_colors()
        section_color = color or C.CYAN
        width = self.terminal_width
        
        print(f"\n{section_color}{'═' * width}{C.END}")
        if title:
            print(f"{section_color}{C.BOLD}{title.center(width)}{C.END}")
            print(f"{section_color}{'═' * width}{C.END}\n")
    
    def display_warning(self, message):
        """Display warning message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.YELLOW}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.YELLOW}{C.BOLD}⚠️  WARNING:{C.END} {C.YELLOW}{message}{C.END}")
        else:
            print(f"{C.YELLOW}{C.BOLD}⚠️  WARNING:{C.END}")
            for line in message_lines:
                print(f"{C.YELLOW}  {line}{C.END}")
        print(f"{C.YELLOW}{'─' * width}{C.END}")
    
    def display_error(self, message):
        """Display error message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.RED}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.RED}{C.BOLD}❌ ERROR:{C.END} {C.RED}{message}{C.END}")
        else:
            print(f"{C.RED}{C.BOLD}❌ ERROR:{C.END}")
            for line in message_lines:
                print(f"{C.RED}  {line}{C.END}")
        print(f"{C.RED}{'─' * width}{C.END}")
    
    def display_success(self, message):
        """Display success message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.GREEN}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.GREEN}{C.BOLD}✅ SUCCESS:{C.END} {C.GREEN}{message}{C.END}")
        else:
            print(f"{C.GREEN}{C.BOLD}✅ SUCCESS:{C.END}")
            for line in message_lines:
                print(f"{C.GREEN}  {line}{C.END}")
        print(f"{C.GREEN}{'─' * width}{C.END}")
    
    def display_info(self, message):
        """Display info message"""
        C = self._get_colors()
        width = self.terminal_width
        
        # Wrap long messages
        wrapper = textwrap.TextWrapper(width=width-4, initial_indent='  ', subsequent_indent='  ')
        message_lines = wrapper.wrap(message)
        
        print(f"{C.CYAN}{'─' * width}{C.END}")
        if len(message_lines) == 1:
            print(f"{C.CYAN}{C.BOLD}💡 INFO:{C.END} {C.CYAN}{message}{C.END}")
        else:
            print(f"{C.CYAN}{C.BOLD}💡 INFO:{C.END}")
            for line in message_lines:
                print(f"{C.CYAN}  {line}{C.END}")
        print(f"{C.CYAN}{'─' * width}{C.END}")
    
    def display_progress_bar(self, iteration, total, prefix='', suffix='', length=50, fill='█'):
        """Display enhanced progress bar"""
        C = self._get_colors()
        
        # Don't display if we're at 0% or 100% and suffix is empty
        if iteration == 0 and not suffix:
            return
            
        percent = ("{0:.1f}").format(100 * (iteration / float(total)))
        filled_length = int(length * iteration // total)
        bar = fill * filled_length + '░' * (length - filled_length)
        
        # Color gradient based on progress
        if iteration < total * 0.33:
            color = C.RED
        elif iteration < total * 0.66:
            color = C.YELLOW
        else:
            color = C.GREEN
        
        # For intermediate updates, use \r to overwrite
        print(f'\r{C.CYAN}{prefix}{C.END} {color}[{bar}]{C.END} {C.BOLD}{percent}%{C.END} {C.DIM}{suffix}{C.END}', end='\r')
        
        # If this is the final iteration, print newline
        if iteration == total:
            print()
            sys.stdout.flush()

def main():
    """Test the banner display with the 8-step startup"""
    banner = StaphScopeBanner()
    
    # Display full banner with startup sequence (8 steps)
    banner.display_startup_sequence()
    banner.display_banner(show_quote=True, show_author=True)
    
    # Test module headers
    banner.display_module_header("MLST Analysis", "Multi-Locus Sequence Typing for S. aureus")
    banner.display_info("Copied 1 files to MLST module")
    banner.display_info("Running MLST analysis with pattern: '*.fna'")
    banner.display_success("MLST analysis completed!")
    
    print()
    banner.display_module_header("spa Typing Analysis", "Staphylococcal Protein A Typing")
    banner.display_info("Copied 1 files to spa module")
    banner.display_success("spa typing analysis completed!")
    
    # Test progress bar
    print("\n🔄 Progress demonstration:")
    for i in range(101):
        banner.display_progress_bar(i, 100, prefix='Analyzing:', suffix='Complete', length=50)
        time.sleep(0.02)
    
    # Display footer
    banner.display_footer(analysis_time="2 minutes, 15 seconds", samples_processed=5)

if __name__ == "__main__":
    main()