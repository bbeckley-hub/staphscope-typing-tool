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
import os
import shutil

class StaphScopeBanner:
    """StaphScope Banner Display with Scientific Quotes"""
    
    def __init__(self):
        self.banner_art = self._get_banner_art()
        self.quotes = self._get_scientific_quotes()
        self.version = "v1.0.0"
        self.author_info = self._get_author_info()
    
    def _get_optimal_width(self, content_type="general"):
        """Get optimal width based on content type"""
        try:
            terminal_width = shutil.get_terminal_size().columns
            
            if content_type == "quote":
                # Quotes look best at moderate width for readability
                return min(max(90, terminal_width - 20), 110)
            elif content_type == "author":
                # Author info doesn't need to be too wide
                return min(max(85, terminal_width - 25), 100)
            elif content_type == "footer":
                # Footer can be a bit wider for the completion message
                return min(max(95, terminal_width - 15), 120)
            else:
                # Default for modules
                return 80
        except:
            return {
                "quote": 100,
                "author": 90, 
                "footer": 100,
                "general": 80
            }[content_type]
    
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
    
    def _create_box(self, width, char='═', corner_left='╔', corner_right='╗'):
        """Create a horizontal box line"""
        return f"{corner_left}{char * (width - 2)}{corner_right}"
    
    def _create_text_line(self, text, width, color=None, align='left', padding=2):
        """Create a text line within the box"""
        C = self._get_colors()
        text_color = color or C.WHITE
        
        # Calculate available width (subtract borders and padding)
        available_width = width - (padding * 2) - 2  # -2 for borders
        
        # Handle text wrapping
        if len(text) > available_width:
            text = textwrap.fill(text, width=available_width)
        
        lines = text.split('\n') if '\n' in text else [text]
        
        formatted_lines = []
        for line in lines:
            if align == 'center':
                formatted_text = line.center(available_width)
            elif align == 'right':
                formatted_text = line.rjust(available_width)
            else:  # left align
                formatted_text = line.ljust(available_width)
            
            formatted_lines.append(f"║{' ' * padding}{text_color}{formatted_text}{C.END}{' ' * padding}║")
        
        return formatted_lines
    
    def display_banner(self, show_quote=True, show_author=True):
        """Display the main StaphScope banner"""
        C = self._get_colors()
        
        print(f"{C.CYAN}{C.BOLD}{self.banner_art}{C.END}")
        print(f"{C.YELLOW}{C.BOLD}         Version: {self.version}{C.END}")
        print()
        
        if show_quote:
            quote = random.choice(self.quotes)
            quote_width = self._get_optimal_width("quote")
            
            # Quote box
            print(f"{C.GREEN}{self._create_box(quote_width, '═', '╔', '╗')}{C.END}")
            
            # Title
            title_lines = self._create_text_line("💡 SCIENTIFIC QUOTE", quote_width, C.CYAN + C.BOLD, 'left')
            for line in title_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(quote_width, '─', '╠', '╣')}{C.END}")
            
            # Quote text
            quote_lines = self._create_text_line(quote['quote'], quote_width, C.WHITE, 'left')
            for line in quote_lines:
                print(line)
            
            # Author
            author_lines = self._create_text_line(f"— {quote['author']}", quote_width, C.YELLOW, 'right')
            for line in author_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(quote_width, '═', '╚', '╝')}{C.END}")
            print()
        
        if show_author:
            author_width = self._get_optimal_width("author")
            
            # Author info box
            print(f"{C.GREEN}{self._create_box(author_width, '═', '╔', '╗')}{C.END}")
            
            # Title
            title_lines = self._create_text_line("👨‍💻 AUTHOR INFORMATION", author_width, C.MAGENTA + C.BOLD, 'left')
            for line in title_lines:
                print(line)
            
            print(f"{C.GREEN}{self._create_box(author_width, '─', '╠', '╣')}{C.END}")
            
            # Author details
            details = [
                f"Name: {self.author_info['name']}",
                f"GitHub: {self.author_info['github']}",
                f"Email: {self.author_info['email']}",
                f"Affiliation: {self.author_info['affiliation']}",
                f"License: {self.author_info['license']}"
            ]
            
            for detail in details:
                detail_lines = self._create_text_line(detail, author_width, C.WHITE, 'left')
                for line in detail_lines:
                    print(line)
            
            print(f"{C.GREEN}{self._create_box(author_width, '═', '╚', '╝')}{C.END}")
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
        
        footer_width = self._get_optimal_width("footer")
        
        print()
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╔', '╗')}{C.END}")
        
        # Title
        title_lines = self._create_text_line("ANALYSIS COMPLETE - STAPHSCOPE", footer_width, C.MAGENTA + C.BOLD, 'center')
        for line in title_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╠', '╣')}{C.END}")
        
        # Analysis details
        if analysis_time:
            time_lines = self._create_text_line(f"⏱️  Analysis Duration: {analysis_time}", footer_width, C.CYAN, 'left')
            for line in time_lines:
                print(line)
        
        if samples_processed > 0:
            samples_lines = self._create_text_line(f"🧫 Samples Processed: {samples_processed}", footer_width, C.GREEN, 'left')
            for line in samples_lines:
                print(line)
        
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        time_lines = self._create_text_line(f"📅 Completion Time: {current_time}", footer_width, C.YELLOW, 'left')
        for line in time_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╠', '╣')}{C.END}")
        
        # Support information
        support_lines = self._create_text_line("TECHNICAL SUPPORT & INQUIRIES:", footer_width, C.WHITE + C.BOLD, 'left')
        for line in support_lines:
            print(line)
        
        author_lines = self._create_text_line(f"Author: {self.author_info['name']}", footer_width, C.CYAN, 'left')
        for line in author_lines:
            print(line)
        
        github_lines = self._create_text_line(f"GitHub: {self.author_info['github']}", footer_width, C.BLUE, 'left')
        for line in github_lines:
            print(line)
        
        email_lines = self._create_text_line(f"Email: {self.author_info['email']}", footer_width, C.GREEN, 'left')
        for line in email_lines:
            print(line)
        
        affiliation_lines = self._create_text_line(f"Affiliation: {self.author_info['affiliation']}", footer_width, C.WHITE, 'left')
        for line in affiliation_lines:
            print(line)
        
        print(f"{C.MAGENTA}{self._create_box(footer_width, '═', '╚', '╝')}{C.END}")
        print()
    
    def display_module_header(self, module_name, description=""):
        """Display header for specific analysis modules"""
        C = self._get_colors()
        
        module_width = self._get_optimal_width("general")
        
        print()
        print(f"{C.BLUE}{self._create_box(module_width, '═', '╔', '╗')}{C.END}")
        
        # Module name
        module_lines = self._create_text_line(f"🧬 {module_name.upper()}", module_width, C.BLUE + C.BOLD, 'left')
        for line in module_lines:
            print(line)
        
        if description:
            print(f"{C.BLUE}{self._create_box(module_width, '─', '╠', '╣')}{C.END}")
            # Description
            desc_lines = self._create_text_line(description, module_width, C.WHITE, 'left')
            for line in desc_lines:
                print(line)
        
        print(f"{C.BLUE}{self._create_box(module_width, '═', '╚', '╝')}{C.END}")
        print()
    
    # ... keep the other methods (display_warning, display_error, etc.) the same ...
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
