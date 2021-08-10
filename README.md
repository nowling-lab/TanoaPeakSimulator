# Tanoa Peak Simulator And Caller

---
Tanoa is the first step of a two step process of peak calling. Tanoa's role in this is to generate a file
of potential peaks to then be filtered through.
---
# Documentation
---
Tutorials and installation guides 
    1. Installing Tanoa (Get it on pip so that they can just use that. No need to clone anything)
    2. Tanoa call peaks
    3. Tanoa generate peaks

---
Quick Start Guide
---
   It is recommended to run Tanoa in a virtual python environment. Information on how to do that is here(https://docs.python.org/3/tutorial/venv.html)

   tanoa_generate_peaks: It's recommended to use the default enhancer length and read length. Needed flags
   are the number of enhancers, the number of reads, the sequence file, and the directory to output.
   An example usage of this program is this:
        
    tanoa_generate_peaks --num-enhancers 3 --num-reads 50 --sequence-file ~/sequence-file.fasta --outputdir ~/outputdirectory

   3 enhancers with 50 reads each generate decent peak-like structures, while also managing a relatively low (compared to higher numbers) run-time

   This will generate 3 regions which will have a peak-like structure each. These are potentially useful in testing for weird peak shapes, as it's generally more simple to use      real data for general usage and testing.
        
   tanoa_call_peaks: This will call peaks from a given file formatted in:
   chromosome, fragment_location, fragment_depth format. An example of this is this: 2L,1632257,1   
   Tanoa compresses this data to use for peak calling, if you want to see that data you may use the --write-depths flag. 
    
   Example usage with and without this flag for call peaks is as such:
            
        tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory
            and
        tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory --write-depths

---

Note: In depth Tanoa notes on the processes of what tanoa does and why can be found here, in the docs folder,
    or in the respective guides for call peak or generate peaks
