The basis of the methodology behind how Tanoa calls peaks is based on the process
of keeping track of the way in which depths change over time. A purely calculus
approach would be the academically ideal for this task in finding where peaks
within a file occur. However, because of the non-standard and relatively
low "resolution" (meaning, the low number of discrete Y steps for most peaks) of
peaks in general, this approach falls short.

To combat this, a new algorithm based on compressing and following how the depths
change was developed. The broad steps are as follows:

1: Compression
    A depth file is a position: depth formatted file in which each position is
    given a number of how many reads were given at that position.

    Examples:
      This sequence of depths:
      2L	5696	3
      2L	5697	3
      2L	5698	3
      2L	5699	3
      2L	5700	3
      2L	5701	3
      Would be compressed down to 5696, 5701: 3. This gives a much more reasonable
      file size (~roughly 100 time smaller) for iteration, and also gives a much
      broader overview of the overall behavior of the depths over time.
      2L 4920 5696 1
      2L 5696 5827 3
      2L 5827 5831 10
      2L 5831 5833 11
      2L 5833 5834 12
      2L 5834 6540 13
      2L 6540 6541 12
      2L 6541 7768 3
      2L 7768 8236 1
      An example of a small peak in a compressed depth file. This also displays
      the low resolution of this peak. Though, in this case a calculus based
      approach would have registered this peak because of how it has two 12's
      surrounding the maximum depth of 13.

2: Rough peak calling:
    The next step is roughly calling the peak by tracking the shift from
    increasing to decreasing. This step registers and calls all peaks by tracking
    each increasing depth and checking the local area of 3 steps(discrete heights in a
    compressed depth file) to the left and right of a point and registering if that
    point is a local maximum, or peak. This step generally has extra called peaks on
    each visual peak and as such is considered noisy, which leads to the next step

3: Cleanup:
    After calling peaks, these peaks are cleaned up. Cleaning up peaks is done by
    finding and locating a background level of depths around a given peak point,
    and combining each and every peak within that range. With the previous
    compressed depth file output as an example, background would be registered
    at 5696, 5827 and 6541, 7768 and all peaks found within those bounds would be
    combined into a singular peak. This maintains accuracy, and gets rid of a lot of
    noise within a peak file.
