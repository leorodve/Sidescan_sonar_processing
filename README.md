# Side Scan Sonar processing and automatic interpreting tool
This is a tool developed as part of my master thesis research. It's divided in two part, the first one includes a basic GUI to allow the user to upload as many sonar files (must be in '.xtf' format) as he desires to process/visualize the data. The second part was written using jupyter notebooks, and it's designed to automatically clasify Side Scan Sonar mosaics into three different facies: ripples, not-ripples and transition zones (tiles that include ripple and non ripple areas). This second part it's also designed to do boulder detection on the mosaics and output a list of the detected boulders, including information of each item such as height, length, confidence level of detection, among many others parameters.

This tool is licensed under the GNU Lesser General Public License (LGPL) v3.0.

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

## Tutorial
You will need to clone this repository to your computer, you can do so by running the following line in your console:

` git clone https://github.com/leorodve/Sidescan_sonar_processing `

You can now run the 'main.py' file, after which you will be greeted with the welcoming window where you need to input some of the navigation and processing parameters or use the default values. You will also need to upload your side scan sonar files in this window.

![Screenshot of welcome window](https://imgur.com/fGjtPzL.png)

If no errors are encountered, you will be able to visualize your data in bird-view.

![Screenshot of welcome window](https://imgur.com/kYEXjL8.png)

And by clicking the 'Waterfall view' button you get to visualize your data in waterfall mode and change some of the display parameters, as well as apply a slant correction function. Once the data has been slant corrected, you can also apply two different types of gain.

![Screenshot of welcome window](https://imgur.com/NYMzG23.png)

![Screenshot of welcome window](https://imgur.com/QpfS4QQ.png)
