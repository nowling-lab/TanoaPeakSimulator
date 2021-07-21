# TanoaPeakSimulator
Basic Usage:

Pre-conditions:

    1: Run a Python Virtual Environment: https://docs.python.org/3/tutorial/venv.html

    2: Run "python3 setup.py install" in the TanoaPeakSimulator directory after cloning the repo with git

Usage:
    tanoa_generate_peaks: It's recommended to use the default enhancer length and read length. Needed flags
    are the number of enhancers, the number of reads, the sequence file, and the directory to output.
    And example usage of this program is this:
    
        tanoa_generate_peaks --num-enhancers 3 --num-reads 50 --sequence-file ~/sequence-file.fasta --outputdir ~/outputdirectory

        3 enhancers with 50 reads each generate decent peak-like structures, while also managing a relatively low (compared to higher numbers) run-time
    
    tanoa_call_peaks: There are no user customizable effects with the call peaks functionality of tanoa. However, there are some diognostic tools that are avaible for the curious user. For an in depth look at 
    the methodology behind how Tanoa calls peak, view the file within docs called Algorithm.md. The compressed depth file that is referenced there may be viewed in its entirety for a given depth file by adding the --write-depths flag. Example usage with and without this flag for call peaks is as such:
        
        tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory

        and

        tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory --write-depths

