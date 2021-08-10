Tanoa Call Peaks
---

Usage
---

The main purpose that Tanoa will be used for is to Call Peaks. Tanoa is built to be the first step in a two 
step process of peak calling, where Tanoa utilizes a very naive approach (detailed here(link)) to peak calling.
This naive approach to peak calling is to simply grab all local maximums within a given file. Tanoa does some
additional processing for peak width and other factors, but Tanoa's main goal is to not miss peaks. 

To use Tanoa to call peaks, you can use the tanoa_call_peaks function with a depth file and a output 
directory. The optional --write-depths flag simply prints the compressed depths (as detailed here(link)) 
for those who are curious.  

For those who do not have a depth file ready, one can be formulated by using the following scripts:
    (link to align_reads_to_fragments.py and calculate_fragment_depth.py)
    (have these stored in the repo somewhere and give tanoa another call which will call them both???)

tanoa_call_peaks:            
            tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory

            and

            tanoa_call_peaks --depth-file ~/fragments_depths.tsv --outputdir ~/outputdirectory --write-depths

