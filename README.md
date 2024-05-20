# TMM-Tool-Acoustic-Modeling
Code implementation of a straightforward, general, efficient and stable TMM proposed by Guochenhao Song, Zhuang Mo, and J. Stuart Bolton

[1] Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton. "A general and stable approach to modeling and coupling multilayered acoustical systems with various types of layers." Journal of Sound and Vibration 567 (2023): 117898.

[2] Song, Guochenhao, Zhuang Mo, and J. Stuart Bolton, "A transfer-matrix-based approach to predicting acoustic properties of a layered system in a general, efficient, and stable Way," SAE Int. J. Adv. & Curr. Prac. in Mobility 6(2):922-934, 2024.



Layered materials are one of the most commonly used acoustical treatments. This tool implemented a method to model and couple layered systems with various layer types (i.e., poro-elastic lavers, solid-elastic layers, stiff panels, and fluid layers) and to stably predict their acoustical properties

This method involves only the topmost layer and its boundary conditions at two interfaces at a time which are further simplified into an equivalent interface. As a result, for a multi-layered system, the proposed method is computationally less expensive.

Moreover, the accuracy of the wave attenuation terms is ensured by decomposing each layer's transfer matrix analytically and reformulating the equation system.
Therefore, this method can produce a stable prediction of acoustical properties over a large frequency and parameter region.
