\name{solar}
\docType{data}
\alias{solar irradiance}
\alias{solar.hs}
\alias{solar.lean}
\alias{beryllium}
\title{Solar Irradiance Proxy Data}
\description{
solar irradiance proxy data.

Hoyt and Schatten (1993) reconstructed solar irradiance (from 1700 through 1997)
using the amplitude of the 11-year solar cycle together with a long term trend estimated
from solar-like stars. They put relatively more weight on the length of the 11-year cycle.

Lean et al. (1995) reconstructed solar irradiance (from 1610 through 2000)
using the amplitude of the 11-year solar cycle and a long term trend estimated
from solar-like stars.

10-Beryllium (10Be) is measured in polar ice from 1424 through 1985.
10-Beryllium (10Be) is produced in the atmosphere
by incoming cosmic ray flux, which in turn is influenced by
the solar activity. The higher the solar activity, the lower
the flux of cosmic radiation entering the earth atmosphere
and therefore the lower the production rate of 10Be. The
short atmospheric lifetime of 10Be of one to two years (Beer
et al. 1994) allows the tracking of solar activity changes and
offers an alternative way to the sunspot based techniques for
the analysis of the amplitude and length of the solar cycle
as well as for low frequency variations.
}
\usage{
data(solar.hs)
data(solar.lean)
data(beryllium)
}
\format{
A list of year and solar (solar irradiance proxy data) for solar.hs and solar.lean
A list of year and be (10-Beryllium) for beryllium
}

\references{
Beer, J., Baumgartner, S., Dittrich-Hannen, B., Hauenstein, J.,
Kubik, P., Lukasczyk, C., Mende, W., Stellmacher, R. and Suter,
M. (1994) Solar variability traced by cosmogenic isotopes. \emph{In:
Pap, J.M., FrQohlich, C., Hudson, H.S., Solanki, S. (Eds.), The
Sun as a Variable Star: Solar and Stellar Irradiance Variations},
Cambridge University Press, Cambridge, 291--300.

Beer, J., Mende, W. and Stellmacher, R. (2000) The role of the sun in
climate forcing. \emph{Quaternary Science Reviews}, \bold{19}, 403--415.

Hoyt, D. V and, Schatten, K. H. (1993) A discussion of plausible
solar irradiance variations, 1700--1992. \emph{Journal of Geophysical
Research}, \bold{98 (A11)}, 18,895--18,906.

Lean, J. L., Beer, J. and Bradley, R. S. (1995) Reconstruction of
solar irradiance since 1610: Implications for climate change.
\emph{Geophysical Research Letters}, \bold{22 (23)}, 3195--3198.

Oh, H-S, Ammann, C. M., Naveau, P., Nychka, D. and Otto-Bliesner, B. L. 
(2003) Multi-resolution time series analysis applied to solar irradiance
and climate reconstructions. 
\emph{Journal of Atmospheric and Solar-Terrestrial Physics}, \bold{65}, 191--201.
}

\examples{
data(solar.hs)
names(solar.hs)
plot(solar.hs$year, solar.hs$solar, type="l")

data(solar.lean)
names(solar.lean)
plot(solar.lean$year, solar.lean$solar, type="l")

data(beryllium)
names(beryllium)
plot(beryllium$year, beryllium$be, type="l")
}
\keyword{datasets}
