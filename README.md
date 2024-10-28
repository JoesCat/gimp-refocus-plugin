# Refocus

- Introduction
- Documentation
- Examples
- News
- Links
- Author

## Introduction

During image processing operations such as scanning and scaling, images tend to get blurry.

The blurred impression of these images is due to the fact that image pixels are averaged with their neighbors. Blurred images don't have sharp boundaries and look as though they have been taken with an unfocussed camera.

Refocus is a plug-in for the [GIMP](http://www.gimp.org/) (the GNU Image Manipulation Program). This plug-in attempts to "refocus" the image, using a technique called FIR Wiener filtering. The traditional technique for sharpening images is to use unsharp masking. Refocus generally produces better results than unsharp masking.

The plug-in comes with a preview that helps you select the best parameters.

The plugin can be found under "Filters > Enhance > Refocus"

## Documentation

Instructions for installing and using the plug-in are included in the distribution in both html and pdf format. An online html version can be found [here](https://refocus.sourceforge.net/doc.html).

## Examples

The following images give you an impression of the effect of this plug-in.

![](img/squirrel.jpg)

Here is the refocussed image of the squirrel

![](img/squirrel-refocussed.jpg)

This Wilber image (unfocussed)

![](img/wilber-unsharp.png)

and optimal refocussing

![](img/wilber-refocussed.png)

## News

- V-0-9-1, 2014 Apr 7 - Additional patches by Ernst Lippe 2004, and added patches by Gentoo, copied from [FreeBSD 2014](https://www.freshports.org/graphics/gimp-refocus-plugin/).
- V-0-9-0, 2003 feb 3, Ernst Lippe - The last published version of Refocus.

## Links

The Refocus home page (2004) was https://refocus.sourceforge.net. The Sourceforge project home page was https://sourceforge.net/projects/refocus. Code was copied and converted from cvs to git format.

## Author

This plug-in was written by Ernst Lippe (2003, 2004).
