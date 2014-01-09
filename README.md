DiscreteSeries
==============

Written by Paul Schellin at the Ohio State University, 2012 to 2013.


DiscreteSeries is a C++ header-only library which encapsulates the relation between a domain and a codomain represented by a series of data points.

Both consistent and inconsistent sampling rates are supported at the moment (and on that note, pardon the number of constructors the class has!).


An interpolation interface is outlined in the FakeInterpolator class, and any class satisfying this interface can be used to interpolate values when a domain value request is made.


The intention is to give this class complete unit support (units as in units of measure, physical quantities) which will most likely be fulfilled through the Boost Units library.


I attempted to make this class as close to the STL container specifications, but because of the nature of encapsulating two domains leads to it being fairly different conceptually from all of the other STL containers (except std::map), I was unable to match the spec entirely.

(More updates on the way!)
