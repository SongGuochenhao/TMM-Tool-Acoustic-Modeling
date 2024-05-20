# TMM-Tool-Acoustic-Modeling
Code implementation of a straightforward, general, efficient and stable TMM proposed by Guochenhao Song, Zhuang Mo, and J. Stuart Bolton




Layered materials are one of the most commonly used acoustical treatments. This tool implemented a method to model and couple layered systems with various layer types (i.e., poro-elastic lavers, solid-elastic layers, stiff panels, and fluid layers) and to stably predict their acoustical properties

This method involves only the topmost layer and its boundary conditions at two interfaces at a time which are further simplified into an equivalent interface. As a result, for a multi-layered system, the proposed method is computationally less expensive.

Moreover, the accuracy of the wave attenuation terms is ensured by decomposing each layer's transfer matrix analytically and reformulating the equation system.
Therefore, this method can produce a stable prediction of acoustical properties over a large frequency and parameter region.
