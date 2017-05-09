## Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
## Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
## Copyright (C) 2011 Alexander Klein <alexander.klein@math.uni-giessen.de>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{b}, @var{a}] =} butter (@var{n}, @var{w})
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} butter (@var{n}, @var{w}, "high")
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} butter (@var{n}, [@var{wl}, @var{wh}])
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} butter (@var{n}, [@var{wl}, @var{wh}], "stop")
## @deftypefnx {Function File} {[@var{z}, @var{p}, @var{g}] =} butter (@dots{})
## @deftypefnx {Function File} {[@var{a}, @var{b}, @var{c}, @var{d}] =} butter (@dots{})
## @deftypefnx {Function File} {[@dots{}] =} butter (@dots{}, "s")
## Generate a Butterworth filter.
## Default is a discrete space (Z) filter.
##
## [b,a] = butter(n, Wc)
##    low pass filter with cutoff pi*Wc radians
##
## [b,a] = butter(n, Wc, 'high')
##    high pass filter with cutoff pi*Wc radians
##
## [b,a] = butter(n, [Wl, Wh])
##    band pass filter with edges pi*Wl and pi*Wh radians
##
## [b,a] = butter(n, [Wl, Wh], 'stop')
##    band reject filter with edges pi*Wl and pi*Wh radians
##
## [z,p,g] = butter(...)
##    return filter as zero-pole-gain rather than coefficients of the
##    numerator and denominator polynomials.
##
## [...] = butter(...,'s')
##     return a Laplace space filter, W can be larger than 1.
##
## [a,b,c,d] = butter(...)
##  return  state-space matrices
##
## References:
##
## Proakis & Manolakis (1992). Digital Signal Processing. New York:
## Macmillan Publishing Company.
## @end deftypefn

function [a, b, c, d] = butter (n, w, varargin)

  if (nargin > 4 || nargin < 2 || nargout > 4 || nargout < 2)
    print_usage ();
  endif

  ## interpret the input parameters
  if (! (isscalar (n) && (n == fix (n)) && (n > 0)))
    error ("butter: filter order N must be a positive integer");
  endif

  stop = false;
  digital = true;
  for i = 1:numel (varargin)
    switch (varargin{i})
      case "s"
        digital = false;
      case "z"
        digital = true;
      case {"high", "stop"}
        stop = true;
      case {"low", "pass"}
        stop = false;
      otherwise
        error ("butter: expected [high|stop] or [s|z]");
    endswitch
  endfor

  if (! ((numel (w) <= 2) && (rows (w) == 1 || columns (w) == 1)))
    error ("butter: frequency must be given as WC or [WL, WH]");
  elseif ((numel (w) == 2) && (w(2) <= w(1)))
    error ("butter: W(1) must be less than W(2)");
  endif

  if (digital && ! all ((w >= 0) & (w <= 1)))
    error ("butter: all elements of W must be in the range [0,1]");
  elseif (! digital && ! all (w >= 0))
    error ("butter: all elements of W must be in the range [0,inf]");
  endif

  ## Prewarp to the band edges to s plane
  if (digital)
    T = 2;       # sampling frequency of 2 Hz
    w = 2 / T * tan (pi * w / T);
  endif

  ## Generate splane poles for the prototype Butterworth filter
  ## source: Kuc
  C = 1;  ## default cutoff frequency
  pole = C * exp (1i * pi * (2 * [1:n] + n - 1) / (2 * n));
  if (mod (n, 2) == 1)
    pole((n + 1) / 2) = -1;  ## pure real value at exp(i*pi)
  endif
  zero = [];
  gain = C^n;

  ## splane frequency transform
  [zero, pole, gain] = sftrans (zero, pole, gain, w, stop);

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
    [zero, pole, gain] = bilinear (zero, pole, gain, T);
  endif

  ## convert to the correct output form
  if (nargout == 2)
    a = real (gain * poly (zero));
    b = real (poly (pole));
  elseif (nargout == 3)
    a = zero;
    b = pole;
    c = gain;
  else
    ## output ss results
    [a, b, c, d] = zp2ss (zero, pole, gain);
  endif

endfunction

%!shared sf, sf2, off_db
%! off_db = 0.5;
%! ## Sampling frequency must be that high to make the low pass filters pass.
%! sf = 6000; sf2 = sf/2;
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];

%!test
%! ## Test low pass order 1 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 1, 50 / sf2 );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 4 ) - damp_db( 5 ), damp_db( 1 : 3 ) ], [ 6 0 0 -3 ], off_db )

%!test
%! ## Test low pass order 4 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 4, 50 / sf2 );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 4 ) - damp_db( 5 ), damp_db( 1 : 3 ) ], [ 24 0 0 -3 ], off_db )

%!test
%! ## Test high pass order 1 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 1, 50 / sf2, "high" );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 2 ) - damp_db( 1 ), damp_db( 3 : end ) ], [ 6 -3 0 0 ], off_db )

%!test
%! ## Test high pass order 4 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 4, 50 / sf2, "high" );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 2 ) - damp_db( 1 ), damp_db( 3 : end ) ], [ 24 -3 0 0 ], off_db )

%% Test input validation
%!error [a, b] = butter ()
%!error [a, b] = butter (1)
%!error [a, b] = butter (1, 2, 3, 4, 5)
%!error [a, b] = butter (.5, .2)
%!error [a, b] = butter (3, .2, "invalid")

%!demo
%! sf = 800; sf2 = sf/2;
%! data=[[1;zeros(sf-1,1)],sinetone(25,sf,1,1),sinetone(50,sf,1,1),sinetone(100,sf,1,1)];
%! [b,a]=butter ( 1, 50 / sf2 );
%! filtered = filter(b,a,data);
%!
%! clf
%! subplot ( columns ( filtered ), 1, 1)
%! plot(filtered(:,1),";Impulse response;")
%! subplot ( columns ( filtered ), 1, 2 )
%! plot(filtered(:,2),";25Hz response;")
%! subplot ( columns ( filtered ), 1, 3 )
%! plot(filtered(:,3),";50Hz response;")
%! subplot ( columns ( filtered ), 1, 4 )
%! plot(filtered(:,4),";100Hz response;")

