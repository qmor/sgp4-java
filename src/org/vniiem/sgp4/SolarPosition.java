package org.vniiem.sgp4;

import java.time.LocalDateTime;

import javax.vecmath.Vector4d;

public class SolarPosition {

Eci FindPosition(final LocalDateTime dt)
{
	 double d =Utils.toJulianDay(Utils.getDoubleUnixTime(dt));
	final double mjd = d - Globals.kEPOCH_JAN1_12H_1900;
	
	final double year = 1900 + mjd / 365.25;
	final double T = (mjd + Delta_ET(year) / Globals.kSECONDS_PER_DAY) / 36525.0;
	final double M = Utils.DegreesToRadians(Utils.Wrap360(358.47583
                + Utils.Wrap360(35999.04975 * T)
                - (0.000150 + 0.0000033 * T) * T * T));
	final double L = Utils.DegreesToRadians(Utils.Wrap360(279.69668
                + Utils.Wrap360(36000.76892 * T)
                + 0.0003025 * T*T));
	final double e = 0.01675104 - (0.0000418 + 0.000000126 * T) * T;
	final double C = Utils.DegreesToRadians((1.919460
                - (0.004789 + 0.000014 * T) * T) * Math.sin(M)
                + (0.020094 - 0.000100 * T) * Math.sin(2 * M)
                + 0.000293 * Math.sin(3 * M));
	final double O = Utils.DegreesToRadians(
            Utils.Wrap360(259.18 - 1934.142 * T));
	final double Lsa = Utils.WrapTwoPI(L + C
            - Utils.DegreesToRadians(0.00569 - 0.00479 * Math.sin(O)));
	final double nu = Utils.WrapTwoPI(M + C);
    double R = 1.0000002 * (1 - e * e) / (1 + e * Math.cos(nu));
    final double eps = Utils.DegreesToRadians(23.452294 - (0.0130125
                + (0.00000164 - 0.000000503 * T) * T) * T + 0.00256 * Math.cos(O));
    R = R * Globals.kAU;

    Vector4d solar_position = new Vector4d(R * Math.cos(Lsa),
            R * Math.sin(Lsa) * Math.cos(eps),
            R * Math.sin(Lsa) * Math.sin(eps),
            R);

    return new Eci(dt, solar_position);
}

double Delta_ET(double year)
{
    return 26.465 + 0.747622 * (year - 1950) + 1.886913
        * Math.sin(Globals.kTWOPI * (year - 1975) / 33);
}

}
