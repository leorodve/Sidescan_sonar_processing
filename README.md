# Side Scan Sonar processing and automatic interpreting tool
This is a tool developed as part of my master thesis research. It's divided in two part, the first one includes a basic GUI to allow the user to upload as many sonar files (must be in '.xtf' format) as he desires to process/visualize the data. The second part was written using jupyter notebooks, and it's designed to automatically clasify Side Scan Sonar mosaics into three different facies: ripples, not-ripples and transition zones (tiles that include ripple and non ripple areas), using amplitude and frequency attributes calculated from the input data. This second part it's also designed to do boulder detection on the mosaics and output a list of the detected boulders, including information of each item such as height, length, confidence level of detection, among many others parameters.

## Libraries requirement
Many of the libraries required already come installed with the Anaconda distribution package, but in case you are not using Anaconda, here are the libraries you need to have installed:

* Numpy 1.26.0
* Matplotlib 3.8.0
* Tkinter 8.6.12
* SciPy  1.11.3
* TQDM  4.66.1
* Seaborn  0.13.0
* Pyproj  3.6.1
* PyXTF  1.3.2
* OpenCV  4.8.1.78
* Rasterio  1.3.8

> [!NOTE]
> The versions of the libraries were the ones used when writing and testing this tool, it may or not work for future versions.

## Examples
### Processing

You will need to clone this repository to your computer, you can do so by running the following line in your console:

` git clone https://github.com/leorodve/Sidescan_sonar_processing `


You can now run the 'main_file.py' file, after which you will be greeted with the welcoming window where you need to input some of the navigation and processing parameters or use the default values. You will also need to upload your side scan sonar files in this window.

![Screenshot of welcome window](https://imgur.com/fGjtPzL.png)

If no errors are encountered, you will be able to visualize your data in bird-view.

![Bird-view](https://imgur.com/kYEXjL8.png)

And by clicking the 'Waterfall view' button you get to visualize your data in waterfall mode and change some of the display parameters, as well as apply a slant correction function. Once the data has been slant corrected, you can also apply two different types of gain.

![Waterfall view](https://imgur.com/NYMzG23.png)

![Waterfall view slant-corrected](https://imgur.com/QpfS4QQ.png)

Here is a comparison of a section slant corrected (a) and the same section with an AGC gain applied (b)

![Waterfall view AGC gain](https://imgur.com/SRj0fuD.png)

### Automatic Interpretation
This section of the code was developed and tested using Google Colab free resources.

**-Facies classification**

The classification file allows the user to classify tiles within the mosaic into ripples, non-ripples and transition tiles, returning also a certainty map for each of the classes. For this the program classify the tiles from the mosaic image, calculates frequency and amplitude attributes and classify them. The final classification results from the combination of the these three classification steps.

Section of a mosaic to classify:

![Section of a mosaic](https://imgur.com/o9D4Gbo.png)

Tiles automatically classified as a ripple:

![Tiles classified as ripples](https://imgur.com/AaqHugg.png)

Certainty map for the ripples class:

![Certainty map for ripples class](https://imgur.com/FISnrzl.png)

**-Object detection**

The object detection file allows the user to identify boulders from a mosaic view. To avoid the invalid detection of pockmarks as boulders, the user must supply the program with a rotation angle so as to make the nadir in the mosaic as vertical as possible, after this the program identify all the nadir present using a semantic segmentation algorithm and separate the data from left and right side of the nadir, flipping one of the two sides so as to make all the shadow patterns the same. Finally, after using an object detection algorithm for three different tiles sizes, the program combines these results using an IoU parameter and yielding a single bounding box list for all the mosaic.

Boulders manually detected:

![Section of a mosaic for detection](https://imgur.com/3tLtLXi.png)

Boulders automatically detected:

![Section of a mosaic for detection](https://imgur.com/DjCAVNC.png)

## Contributing

All kinds of contributions are welcomed.

## License

All source code is licensed under a MIT license. 
