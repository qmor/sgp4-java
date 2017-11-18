package org.vniiem.sgp4;

import java.time.LocalDateTime;

public class OrbitalElements {
    double mean_anomoly_;
    double ascending_node_;
    double argument_perigee_;
    double eccentricity_;
    double inclination_;
    double mean_motion_;
    double bstar_;
    double recovered_semi_major_axis_;
    double recovered_mean_motion_;
    double perigee_;
    double period_;
    LocalDateTime epoch_;
    
    public LocalDateTime Epoch()
    {
    	return epoch_;
    }
    public double MeanAnomoly()
    {
    	return mean_anomoly_;
    }
    public double BStar()
    {
    	return bstar_;
    }
    public double ArgumentPerigee()
    {
    	return argument_perigee_;
    }
    
    public double AscendingNode()
    {
    	return ascending_node_;
    }
    public double Perigee()
    {
    	return perigee_;
    }
    public double Period()
    {
    	return period_;
    }
    double RecoveredMeanMotion()
    {
    	return recovered_mean_motion_;
    }
    double RecoveredSemiMajorAxis()
    {
    	return recovered_semi_major_axis_;
    }
    double Inclination()
    {
    	return inclination_;
    }
    double Eccentricity()
    {
    	return eccentricity_;
    }
    double MeanMotion()
    {
    	return mean_motion_;
    }
	public OrbitalElements(Tle tle)
	{
	    /*
	     * extract and format tle data
	     */
	    mean_anomoly_ = tle.MeanAnomaly(false);
	    ascending_node_ = tle.RightAscendingNode(false);
	    argument_perigee_ = tle.ArgumentPerigee(false);
	    eccentricity_ = tle.Eccentricity();
	    inclination_ = tle.Inclination(false);
	    mean_motion_ = tle.MeanMotion() * Globals.kTWOPI / Globals.kMINUTES_PER_DAY;
	    bstar_ = tle.BStar();
	    epoch_ = tle.Epoch();

	    /*
	     * recover original mean motion (xnodp) and semimajor axis (aodp)
	     * from input elements
	     */
	    final double a1 = Math.pow(Globals.kXKE / MeanMotion(), Globals.kTWOTHIRD);
	    final double cosio = Math.cos(Inclination());
	    final double theta2 = cosio * cosio;
	    final double x3thm1 = 3.0 * theta2 - 1.0;
	    final double eosq = Eccentricity() * Eccentricity();
	    final double betao2 = 1.0 - eosq;
	    final double betao = Math.sqrt(betao2);
	    final double temp = (1.5 * Globals.kCK2) * x3thm1 / (betao * betao2);
	    final double del1 = temp / (a1 * a1);
	    final double a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)));
	    final double del0 = temp / (a0 * a0);

	    recovered_mean_motion_ = MeanMotion() / (1.0 + del0);
	    /*
	     * alternative way to calculate
	     * doesnt affect final results
	     * recovered_semi_major_axis_ = pow(XKE / RecoveredMeanMotion(), TWOTHIRD);
	     */
	    recovered_semi_major_axis_ = a0 / (1.0 - del0);

	    /*
	     * find perigee and period
	     */
	    perigee_ = (RecoveredSemiMajorAxis() * (1.0 - Eccentricity()) - Globals.kAE) * Globals.kXKMPER;
	    period_ = Globals.kTWOPI / RecoveredMeanMotion();
	}

}
