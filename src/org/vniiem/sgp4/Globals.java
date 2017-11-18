package org.vniiem.sgp4;



public class Globals {
	public static final double kAE = 1.0;
	public static final double kQ0 = 120.0;
	public static final double kS0 = 78.0;
	public static final double kMU = 398600.8;
	public static final double kXKMPER = 6378.135;
	public static final double kXJ2 = 1.082616e-3;
	public static final double kXJ3 = -2.53881e-6;
	public static final double kXJ4 = -1.65597e-6;

	/*
	 * alternative XKE
	 * affects final results
	 * aiaa-2006-6573
	 * public final double kXKE = 60.0 / sqrt(kXKMPER * kXKMPER * kXKMPER / kMU);
	 * dundee
	 * public final double kXKE = 7.43669161331734132e-2;
	 */
	public static final double kXKE = 60.0 / Math.sqrt(kXKMPER * kXKMPER * kXKMPER / kMU);
	public static final double kCK2 = 0.5 * kXJ2 * kAE * kAE;
	public static final double kCK4 = -0.375 * kXJ4 * kAE * kAE * kAE * kAE;

	/*
	 * alternative QOMS2T
	 * affects final results
	 * aiaa-2006-6573
	 * #define QOMS2T   (pow(((Q0 - S0) / XKMPER), 4.0))
	 * dundee
	 * #define QOMS2T   (1.880279159015270643865e-9)
	 */
	public static final double kQOMS2T = Math.pow(((kQ0 - kS0) / kXKMPER), 4.0);

	public static final double kS = kAE * (1.0 + kS0 / kXKMPER);
	public final static double kPI = 3.14159265358979323846264338327950288419716939937510582;
	public  final static double kTWOPI = 2.0 * kPI;
	public static final double kTWOTHIRD = 2.0 / 3.0;
	public static final double kTHDT = 4.37526908801129966e-3;
	/*
	 * earth flattening
	 */
	public static final double kF = 1.0 / 298.26;
	/*
	 * earth rotation per sideral day
	 */
	public static final double kOMEGA_E = 1.00273790934;
	public static final double kAU = 1.49597870691e8;

	public static final double kSECONDS_PER_DAY = 86400.0;
	public static final double kMINUTES_PER_DAY = 1440.0;
	public static final double kHOURS_PER_DAY = 24.0;

	// Jan 1.0 1900 = Jan 1 1900 00h UTC
	public static final double kEPOCH_JAN1_00H_1900 = 2415019.5;

	// Jan 1.5 1900 = Jan 1 1900 12h UTC
	public static final double kEPOCH_JAN1_12H_1900 = 2415020.0;

	// Jan 1.5 2000 = Jan 1 2000 12h UTC
	public static final double kEPOCH_JAN1_12H_2000 = 2451545.0;
}
