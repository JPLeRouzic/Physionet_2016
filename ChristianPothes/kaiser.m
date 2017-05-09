## Copyright (C) 1995, 1996, 1997 Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
## Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
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
## @deftypefn  {Function File} {} kaiser (@var{m})
## @deftypefnx {Function File} {} kaiser (@var{m}, @var{beta})
##
## Return the filter coefficients of a Kaiser window of length @var{m}.  The
## Fourier transform of the window has a stop-band attenuation that is derived
## from the parameter @var{beta}.
##
## For the definition of the Kaiser window, see A. V. Oppenheim &
## R. W. Schafer, "Discrete-Time Signal Processing".
##
## The continuous version of width m centered about x=0 is:
##
## @example
## @group
##         besseli(0, beta * sqrt(1-(2*x/m).^2))
## k(x) =  -------------------------------------,  m/2 <= x <= m/2
##                besseli(0, beta)
## @end group
## @end example
##
## @seealso{kaiserord}
## @end deftypefn

function w = kaiser (m, beta = 0.5)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("kaiser: M must be a positive integer");
  elseif (! (isscalar (beta) && isreal (beta)))
    error ("kaiser: BETA must be a real scalar");
  endif

  if (m == 1)
    w = 1;
  else
    N = m - 1;
    k = (0 : N)';
    k = 2 * beta / N * sqrt (k .* (N - k));
    w = besseli (0, k) / besseli (0, beta);
  endif

endfunction

%!demo
%! % use demo("kaiserord");

%!assert (kaiser (1), 1)

%% Test input validation
%!error kaiser ()
%!error kaiser (0.5)
%!error kaiser (-1)
%!error kaiser (ones (1, 4))
%!error kaiser (1, 2, 3)
