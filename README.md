# Refocus

- Introduction
- Documentation
- Examples
- Build and Install
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

## Build and Install

As there are many files in the src directory, it will be easier to build
Refocus using the autoconf tools. You can use gimptool-2.0 to install or
use make install. To build Refocus from Git master requires 2 preparatory
steps:

First, you need to create the `./configure` script if you do not have it yet
```sh
autoreconf -i  (or use 'autoreconf --install --force' for more modern setups)
automake --foreign -Wall
```

Second, you then use the usual steps to compile the Refocus plugin.
Various operating systems and setups will need ./configure options set.
The INSTALLATION file has detailed info for `configure' options.
Example install steps shown below are for Linux:
```sh
./configure --prefix=/usr
make
make check
sudo make install
```
and to uninstall, use `make uninstall`.

If you want to use `gimptool-2.0`, you can select to run either
`make install-bin` for a single local user, or if you have administration
root access, you can run `sudo make install-admin-bin`.

If Gimp is already running in your system, exit and restart Gimp for the
new Refocus plug-in to be detected in menu "Filters > Enhance > Refocus".

## News

- V-0-9-1, 2014 Apr 7 - Additional patches by Ernst Lippe 2004, and added patches by Gentoo, copied from [FreeBSD 2014](https://www.freshports.org/graphics/gimp-refocus-plugin/).
- V-0-9-0, 2003 feb 3, Ernst Lippe - The last published version of Refocus.

## Links

The Refocus home page (2004) was https://refocus.sourceforge.net. The Sourceforge project home page was https://sourceforge.net/projects/refocus. Code was copied and converted from cvs to git format.

## Author

This plug-in was written by Ernst Lippe (2003, 2004).
