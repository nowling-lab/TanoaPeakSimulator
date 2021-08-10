# Tanoa
---

Tanoa can be used to generate peak-like structures. This feature currently is not being actively developed
    and is still very experimental. Currently, peaks will be generated using randomly sampled reads from an
    given randomly gathered enhancer region. It is recommended to use the current defaults, but the peak shape 
    can be changed drastically by either extending the peak region or increasing the amount of reads per 
    enhancer region. 

Examples of the currently supported generate peak commands are here:

            tanoa_generate_peaks --num-enhancers 3 --num-reads 50 --sequence-file ~/sequence-file.fasta --outputdir ~/outputdirectory

            3 enhancers with 50 reads each generate decent peak-like structures, while also managing a relatively low (compared to higher numbers) run-time

If you would like to change the values to see if you can find a configuration that works for your usage, here
is a table of what each command does:

| Command | Description |
| ----------- | ----------- |
| --enhancer-length | This setting is defaulted at 1000. This is a region in which a num-reads value of reads will be sampled from |
| --num-enhancers | This flag chooses how many separate enhancer regions you would like to sample from |
| --read-length | This is the length of each individual read that takes place in the enhancer regions. This defaults at 100, and will drastically change the shape of the enhancer region. |
| --num-reads | The number of those reads that will be taken. 50 is a good starting value |
| --sequence-file | A .fasta file with nucleotide data from DNA sequencing. This is where the enhancer regions will be chosen from |
| --outputdir | This is the directory files will be written to from the output of generating peaks | 

