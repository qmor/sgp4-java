package org.vniiem.sgp4;

import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.temporal.ChronoUnit;

import javax.vecmath.Vector3d;
import javax.vecmath.Vector4d;

class CommonConstants
{
    double cosio;
    double sinio;
    double eta;
    double t2cof;
    double a3ovk2;
    double x1mth2;
    double x3thm1;
    double x7thm1;
    double aycof;
    double xlcof;
    double xnodcf;
    double c1;
    double c4;
    double omgdot; // secular rate of omega (radians/sec)
    double xnodot; // secular rate of xnode (radians/sec)
    double xmdot;  // secular rate of xmo   (radians/sec)
};

class NearSpaceConstants
{
    double c5;
    double omgcof;
    double xmcof;
    double delmo;
    double sinmo;
    double d2;
    double d3;
    double d4;
    double t3cof;
    double t4cof;
    double t5cof;
};

class DeepSpaceConstants
{
    double gsto;
    double zmol;
    double zmos;

    /*
     * lunar / solar constants for epoch
     * applied during DeepSpaceSecular()
     */
    double sse;
    double ssi;
    double ssl;
    double ssg;
    double ssh;
    /*
     * lunar / solar constants
     * used during DeepSpaceCalculateLunarSolarTerms()
     */
    double se2;
    double si2;
    double sl2;
    double sgh2;
    double sh2;
    double se3;
    double si3;
    double sl3;
    double sgh3;
    double sh3;
    double sl4;
    double sgh4;
    double ee2;
    double e3;
    double xi2;
    double xi3;
    double xl2;
    double xl3;
    double xl4;
    double xgh2;
    double xgh3;
    double xgh4;
    double xh2;
    double xh3;
    /*
     * used during DeepSpaceCalcDotTerms()
     */
    double d2201;
    double d2211;
    double d3210;
    double d3222;
    double d4410;
    double d4422;
    double d5220;
    double d5232;
    double d5421;
    double d5433;
    double del1;
    double del2;
    double del3;
    /*
     * whether the deep space orbit is
     * geopotential resonance for 12 hour orbits
     */
    boolean resonance_flag;
    /*
     * whether the deep space orbit is
     * 24h synchronous resonance
     */
    boolean synchronous_flag;
};

class IntegratorValues
{
   public double xndot;
   public double xnddt;
   public double xldot;
};

class IntegratorConstants
{
    /*
     * integrator constants
     */
    double xfact;
    double xlamo;

    /*
     * integrator values for epoch
     */
    IntegratorValues values_0= new IntegratorValues();
};

class IntegratorParams
{
    /*
     * integrator values
     */
    double xli;
    double xni;
    double atime;
    /*
     * itegrator values for current d_atime_
     */
    IntegratorValues values_t;
};

public class SGP4{
	
	
	/*
	 * Copyright 2013 Daniel Warner <contact@danrw.com>
	 *
	 * Licensed under the Apache License, Version 2.0 (the "License");
	 * you may not use this file except in compliance with the License.
	 * You may obtain a copy of the License at
	 *
	 * http://www.apache.org/licenses/LICENSE-2.0
	 *
	 * Unless required by applicable law or agreed to in writing, software
	 * distributed under the License is distributed on an "AS IS" BASIS,
	 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	 * See the License for the specific language governing permissions and
	 * limitations under the License.
	 */
    static final double ZNS = 1.19459E-5;
    static final double C1SS = 2.9864797E-6;
    static final double ZES = 0.01675;
    static final double ZNL = 1.5835218E-4;
    static final double C1L = 4.7968065E-7;
    static final double ZEL = 0.05490;
    static final double ZCOSIS = 0.91744867;
    static final double ZSINI = 0.39785416;
    static final double ZSINGS = -0.98088458;
    static final double ZCOSGS = 0.1945905;
    static final double Q22 = 1.7891679E-6;
    static final double Q31 = 2.1460748E-6;
    static final double Q33 = 2.2123015E-7;
    static final double ROOT22 = 1.7891679E-6;
    static final double ROOT32 = 3.7393792E-7;
    static final double ROOT44 = 7.3636953E-9;
    static final double ROOT52 = 1.1428639E-7;
    static final double ROOT54 = 2.1765803E-9;
	    
	    /*
	     * the constants used
	     */
	    CommonConstants common_consts_;
	    NearSpaceConstants nearspace_consts_;
	    DeepSpaceConstants deepspace_consts_;
	    IntegratorConstants integrator_consts_;
	    IntegratorParams integrator_params_;

	    /*
	     * the orbit data
	     */
	    OrbitalElements elements_;

	    /*
	     * flags
	     */
	    boolean use_simple_model_;
	    boolean use_deep_space_;

	    static final CommonConstants Empty_CommonConstants = new CommonConstants();
	    static final NearSpaceConstants Empty_NearSpaceConstants = new NearSpaceConstants();
	    static final DeepSpaceConstants Empty_DeepSpaceConstants = new DeepSpaceConstants();
	    static final IntegratorConstants Empty_IntegratorConstants = new IntegratorConstants();
	    static final IntegratorParams Empty_IntegratorParams = new IntegratorParams();



	
	
	/*
	 * Copyright 2013 Daniel Warner <contact@danrw.com>
	 *
	 * Licensed under the Apache License, Version 2.0 (the "License");
	 * you may not use this file except in compliance with the License.
	 * You may obtain a copy of the License at
	 *
	 * http://www.apache.org/licenses/LICENSE-2.0
	 *
	 * Unless required by applicable law or agreed to in writing, software
	 * distributed under the License is distributed on an "AS IS" BASIS,
	 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	 * See the License for the specific language governing permissions and
	 * limitations under the License.
	 */





	//final SGP4::Commonfinalants SGP4::Empty_Commonfinalants = SGP4::Commonfinalants();
	//final SGP4::NearSpacefinalants SGP4::Empty_NearSpacefinalants = SGP4::NearSpacefinalants();
	//final SGP4::DeepSpacefinalants SGP4::Empty_DeepSpacefinalants = SGP4::DeepSpacefinalants();
	//final SGP4::Integratorfinalants SGP4::Empty_Integratorfinalants = SGP4::Integratorfinalants();
	//final SGP4::IntegratorParams SGP4::Empty_IntegratorParams = SGP4::IntegratorParams();

	public void SetTle(Tle tle) throws SatelliteException
	{
	    /*
	     * extract and format tle data
	     */
	    elements_ = new OrbitalElements(tle);

	    Initialise();
	}

	void Initialise() throws SatelliteException
	{
	    /*
	     * reset all finalants etc
	     */
	    Reset();

	    /*
	     * error checks
	     */
	    if (elements_.Eccentricity() < 0.0 || elements_.Eccentricity() > 0.999)
	    {
	        throw new SatelliteException("Eccentricity out of range");
	    }

	    if (elements_.Inclination() < 0.0 || elements_.Inclination() > Globals.kPI)
	    {
	        throw new SatelliteException("Inclination out of range");
	    }

	    common_consts_.cosio = Math.cos(elements_.Inclination());
	    common_consts_.sinio = Math.sin(elements_.Inclination());
	    final double theta2 = common_consts_.cosio * common_consts_.cosio;
	    common_consts_.x3thm1 = 3.0 * theta2 - 1.0;
	    final double eosq = elements_.Eccentricity() * elements_.Eccentricity();
	    final double betao2 = 1.0 - eosq;
	    final double betao = Math.sqrt(betao2);

	    if (elements_.Period() >= 225.0)
	    {
	        use_deep_space_ = true;
	    }
	    else
	    {
	        use_deep_space_ = false;
	        use_simple_model_ = false;
	        /*
	         * for perigee less than 220 kilometers, the simple_model flag is set
	         * and the equations are truncated to linear variation in sqrt a and
	         * quadratic variation in mean anomly. also, the c3 term, the
	         * delta omega term and the delta m term are dropped
	         */
	        if (elements_.Perigee() < 220.0)
	        {
	            use_simple_model_ = true;
	        }
	    }

	    /*
	     * for perigee below 156km, the values of
	     * s4 and qoms2t are altered
	     */
	    double s4 = Globals.kS;
	    double qoms24 = Globals.kQOMS2T;
	    if (elements_.Perigee() < 156.0)
	    {
	        s4 = elements_.Perigee() - 78.0;
	        if (elements_.Perigee() < 98.0) 
	        {
	            s4 = 20.0;
	        }
	        qoms24 = Math.pow((120.0 - s4) * Globals.kAE / Globals.kXKMPER, 4.0);
	        s4 = s4 / Globals.kXKMPER + Globals.kAE;
	    }

	    /*
	     * generate finalants
	     */
	    final double pinvsq = 1.0
	        / (elements_.RecoveredSemiMajorAxis()
	                * elements_.RecoveredSemiMajorAxis()
	                * betao2 * betao2);
	    final double tsi = 1.0 / (elements_.RecoveredSemiMajorAxis() - s4);
	    common_consts_.eta = elements_.RecoveredSemiMajorAxis()
	        * elements_.Eccentricity() * tsi;
	    final double etasq = common_consts_.eta * common_consts_.eta;
	    final double eeta = elements_.Eccentricity() * common_consts_.eta;
	    final double psisq = Math.abs(1.0 - etasq);
	    final double coef = qoms24 * Math.pow(tsi, 4.0);
	    final double coef1 = coef / Math.pow(psisq, 3.5);
	    final double c2 = coef1 * elements_.RecoveredMeanMotion()
	        * (elements_.RecoveredSemiMajorAxis()
	        * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
	        + 0.75 * Globals.kCK2 * tsi / psisq * common_consts_.x3thm1
	        * (8.0 + 3.0 * etasq * (8.0 + etasq)));
	    common_consts_.c1 = elements_.BStar() * c2;
	    common_consts_.a3ovk2 = -Globals.kXJ3 / Globals.kCK2 * Globals.kAE * Globals.kAE * Globals.kAE;
	    common_consts_.x1mth2 = 1.0 - theta2;
	    common_consts_.c4 = 2.0 * elements_.RecoveredMeanMotion()
	        * coef1 * elements_.RecoveredSemiMajorAxis() * betao2
	        * (common_consts_.eta * (2.0 + 0.5 * etasq) + elements_.Eccentricity()
	        * (0.5 + 2.0 * etasq)
	        - 2.0 * Globals.kCK2 * tsi / (elements_.RecoveredSemiMajorAxis() * psisq)
	        * (-3.0 * common_consts_.x3thm1 * (1.0 - 2.0 * eeta + etasq
	        * (1.5 - 0.5 * eeta))
	        + 0.75 * common_consts_.x1mth2 * (2.0 * etasq - eeta *
	            (1.0 + etasq)) * Math.cos(2.0 * elements_.ArgumentPerigee())));
	    final double theta4 = theta2 * theta2;
	    final double temp1 = 3.0 * Globals.kCK2 * pinvsq * elements_.RecoveredMeanMotion();
	    final double temp2 = temp1 * Globals.kCK2 * pinvsq;
	    final double temp3 = 1.25 * Globals.kCK4 * pinvsq * pinvsq * elements_.RecoveredMeanMotion();
	    common_consts_.xmdot = elements_.RecoveredMeanMotion() + 0.5 * temp1 * betao *
	            common_consts_.x3thm1 + 0.0625 * temp2 * betao *
	            (13.0 - 78.0 * theta2 + 137.0 * theta4);
	    final double x1m5th = 1.0 - 5.0 * theta2;
	    common_consts_.omgdot = -0.5 * temp1 * x1m5th +
	            0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
	            temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
	    final double xhdot1 = -temp1 * common_consts_.cosio;
	    common_consts_.xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 *
	            (3.0 - 7.0 * theta2)) * common_consts_.cosio;
	    common_consts_.xnodcf = 3.5 * betao2 * xhdot1 * common_consts_.c1;
	    common_consts_.t2cof = 1.5 * common_consts_.c1;

	    if (Math.abs(common_consts_.cosio + 1.0) > 1.5e-12)
	    {
	        common_consts_.xlcof = 0.125 * common_consts_.a3ovk2 * common_consts_.sinio * (3.0 + 5.0 * common_consts_.cosio) / (1.0 + common_consts_.cosio);
	    }
	    else
	    {
	        common_consts_.xlcof = 0.125 * common_consts_.a3ovk2 * common_consts_.sinio * (3.0 + 5.0 * common_consts_.cosio) / 1.5e-12;
	    }

	    common_consts_.aycof = 0.25 * common_consts_.a3ovk2 * common_consts_.sinio;
	    common_consts_.x7thm1 = 7.0 * theta2 - 1.0;

	    if (use_deep_space_)
	    {
	        deepspace_consts_.gsto = Utils.ToGreenwichSiderealTime(elements_.Epoch());

	        DeepSpaceInitialise(eosq, common_consts_.sinio, common_consts_.cosio, betao,
	                theta2, betao2,
	                common_consts_.xmdot, common_consts_.omgdot, common_consts_.xnodot);
	    }
	    else
	    {
	        double c3 = 0.0;
	        if (elements_.Eccentricity() > 1.0e-4)
	        {
	            c3 = coef * tsi * common_consts_.a3ovk2 * elements_.RecoveredMeanMotion() * Globals.kAE *
	                    common_consts_.sinio / elements_.Eccentricity();
	        }

	        nearspace_consts_.c5 = 2.0 * coef1 * elements_.RecoveredSemiMajorAxis() * betao2 * (1.0 + 2.75 *
	                (etasq + eeta) + eeta * etasq);
	        nearspace_consts_.omgcof = elements_.BStar() * c3 * Math.cos(elements_.ArgumentPerigee());

	        nearspace_consts_.xmcof = 0.0;
	        if (elements_.Eccentricity() > 1.0e-4)
	        {
	            nearspace_consts_.xmcof = -Globals.kTWOTHIRD * coef * elements_.BStar() * Globals.kAE / eeta;
	        }

	        nearspace_consts_.delmo = Math.pow(1.0 + common_consts_.eta * (Math.cos(elements_.MeanAnomoly())), 3.0);
	        nearspace_consts_.sinmo = Math.sin(elements_.MeanAnomoly());

	        if (!use_simple_model_)
	        {
	            final double c1sq = common_consts_.c1 * common_consts_.c1;
	            nearspace_consts_.d2 = 4.0 * elements_.RecoveredSemiMajorAxis() * tsi * c1sq;
	            final double temp = nearspace_consts_.d2 * tsi * common_consts_.c1 / 3.0;
	            nearspace_consts_.d3 = (17.0 * elements_.RecoveredSemiMajorAxis() + s4) * temp;
	            nearspace_consts_.d4 = 0.5 * temp * elements_.RecoveredSemiMajorAxis() *
	                    tsi * (221.0 * elements_.RecoveredSemiMajorAxis() + 31.0 * s4) * common_consts_.c1;
	            nearspace_consts_.t3cof = nearspace_consts_.d2 + 2.0 * c1sq;
	            nearspace_consts_.t4cof = 0.25 * (3.0 * nearspace_consts_.d3 + common_consts_.c1 *
	                    (12.0 * nearspace_consts_.d2 + 10.0 * c1sq));
	            nearspace_consts_.t5cof = 0.2 * (3.0 * nearspace_consts_.d4 + 12.0 * common_consts_.c1 *
	                    nearspace_consts_.d3 + 6.0 * nearspace_consts_.d2 * nearspace_consts_.d2 + 15.0 *
	                    c1sq * (2.0 * nearspace_consts_.d2 + c1sq));
	        }
	    }
	}

	public Eci FindPosition(LocalDateTime dt) throws SatelliteException, DecayedException
	{
		LocalDateTime restime = Utils.minusdt(dt, elements_.Epoch());
	    return FindPosition( Utils.TotalMinutes(restime));
	}

	Eci FindPosition(double tsince) throws SatelliteException, DecayedException
	{
	    if (use_deep_space_)
	    {
	        return FindPositionSDP4(tsince);
	    }
	    else
	    {
	        return FindPositionSGP4(tsince);
	    }
	}

	Eci FindPositionSDP4(double tsince) throws SatelliteException, DecayedException
	{
	    /*
	     * the final values
	     */
	    double[] e = new double[]{0.0};
	    double a;
	    double omega;
	    double xl;
	    double[] xnode = new double[]{0.0};
	    double[] xincl = new double[]{0.0};

	    /*
	     * update for secular gravity and atmospheric drag
	     */
	    double[] xmdf = new double[] {elements_.MeanAnomoly()
	        + common_consts_.xmdot * tsince};
	    double[] omgadf = new double[] {elements_.ArgumentPerigee()
	        + common_consts_.omgdot * tsince};
	    final double[] xnoddf = new double[] {elements_.AscendingNode()
	        + common_consts_.xnodot * tsince};

	    final double tsq = tsince * tsince;
	    xnode[0] = xnoddf[0] + common_consts_.xnodcf * tsq;
	    double tempa = 1.0 - common_consts_.c1 * tsince;
	    double tempe = elements_.BStar() * common_consts_.c4 * tsince;
	    double templ = common_consts_.t2cof * tsq;

	    double[] xn = new double[] { elements_.RecoveredMeanMotion()};
	    e[0] = elements_.Eccentricity();
	    xincl[0] = elements_.Inclination();

	    DeepSpaceSecular(tsince, xmdf, omgadf, xnode, e, xincl, xn);

	    if (xn[0] <= 0.0)
	    {
	        throw new SatelliteException("Error: (xn <= 0.0)");
	    }

	    a = Math.pow(Globals.kXKE / xn[0], Globals.kTWOTHIRD) * tempa * tempa;
	    e[0] -= tempe;
	    double[] xmam = new double[] {xmdf[0] + elements_.RecoveredMeanMotion() * templ};

	    DeepSpacePeriodics(tsince, e, xincl, omgadf, xnode, xmam);

	    /*
	     * keeping xincl positive important unless you need to display xincl
	     * and dislike negative inclinations
	     */
	    if (xincl[0] < 0.0)
	    {
	        xincl[0] = -xincl[0];
	        xnode[0] += Globals.kPI;
	        omgadf[0] -= Globals.kPI;
	    }

	    xl = xmam[0] + omgadf[0] + xnode[0];
	    omega = omgadf[0];

	    /*
	     * fix tolerance for error recognition
	     */
	    if (e[0] <= -0.001)
	    {
	        throw new SatelliteException("Error: (e <= -0.001)");
	    }
	    else if (e[0] < 1.0e-6)
	    {
	        e[0] = 1.0e-6;
	    }
	    else if (e[0] > (1.0 - 1.0e-6))
	    {
	        e[0] = 1.0 - 1.0e-6;
	    }

	    /*
	     * re-compute the perturbed values
	     */
	    final double perturbed_sinio = Math.sin(xincl[0]);
	    final double perturbed_cosio = Math.cos(xincl[0]);

	    final double perturbed_theta2 = perturbed_cosio * perturbed_cosio;

	    final double perturbed_x3thm1 = 3.0 * perturbed_theta2 - 1.0;
	    final double perturbed_x1mth2 = 1.0 - perturbed_theta2;
	    final double perturbed_x7thm1 = 7.0 * perturbed_theta2 - 1.0;

	    double perturbed_xlcof;
	    if (Math.abs(perturbed_cosio + 1.0) > 1.5e-12)
	    {
	        perturbed_xlcof = 0.125 * common_consts_.a3ovk2 * perturbed_sinio
	            * (3.0 + 5.0 * perturbed_cosio) / (1.0 + perturbed_cosio);
	    }
	    else
	    {
	        perturbed_xlcof = 0.125 * common_consts_.a3ovk2 * perturbed_sinio
	            * (3.0 + 5.0 * perturbed_cosio) / 1.5e-12;
	    }

	    final double perturbed_aycof = 0.25 * common_consts_.a3ovk2
	        * perturbed_sinio;

	    /*
	     * using calculated values, find position and velocity
	     */
	    return CalculateFinalPositionVelocity(tsince, e[0],
	            a, omega, xl, xnode[0],
	            xincl[0], perturbed_xlcof, perturbed_aycof,
	            perturbed_x3thm1, perturbed_x1mth2, perturbed_x7thm1,
	            perturbed_cosio, perturbed_sinio);

	}

	Eci FindPositionSGP4(double tsince) throws SatelliteException, DecayedException
	{
	    /*
	     * the final values
	     */
	    double e;
	    double a;
	    double omega;
	    double xl;
	    double xnode;
	    double xincl;

	    /*
	     * update for secular gravity and atmospheric drag
	     */
	    final double xmdf = elements_.MeanAnomoly()
	        + common_consts_.xmdot * tsince;
	    final double omgadf = elements_.ArgumentPerigee()
	        + common_consts_.omgdot * tsince;
	    final double xnoddf = elements_.AscendingNode()
	        + common_consts_.xnodot * tsince;

	    final double tsq = tsince * tsince;
	    xnode = xnoddf + common_consts_.xnodcf * tsq;
	    double tempa = 1.0 - common_consts_.c1 * tsince;
	    double tempe = elements_.BStar() * common_consts_.c4 * tsince;
	    double templ = common_consts_.t2cof * tsq;

	    xincl = elements_.Inclination();
	    omega = omgadf;
	    double xmp = xmdf;

	    if (!use_simple_model_)
	    {
	        final double delomg = nearspace_consts_.omgcof * tsince;
	        final double delm = nearspace_consts_.xmcof
	            * (Math.pow(1.0 + common_consts_.eta * Math.cos(xmdf), 3.0)
	                    * - nearspace_consts_.delmo);
	        final double temp = delomg + delm;

	        xmp += temp;
	        omega -= temp;

	        final double tcube = tsq * tsince;
	        final double tfour = tsince * tcube;

	        tempa = tempa - nearspace_consts_.d2 * tsq - nearspace_consts_.d3
	            * tcube - nearspace_consts_.d4 * tfour;
	        tempe += elements_.BStar() * nearspace_consts_.c5
	            * (Math.sin(xmp) - nearspace_consts_.sinmo);
	        templ += nearspace_consts_.t3cof * tcube + tfour
	            * (nearspace_consts_.t4cof + tsince * nearspace_consts_.t5cof);
	    }

	    a = elements_.RecoveredSemiMajorAxis() * tempa * tempa;
	    e = elements_.Eccentricity() - tempe;
	    xl = xmp + omega + xnode + elements_.RecoveredMeanMotion() * templ;

	    /*
	     * fix tolerance for error recognition
	     */
	    if (e <= -0.001)
	    {
	        throw new SatelliteException("Error: (e <= -0.001)");
	    }
	    else if (e < 1.0e-6)
	    {
	        e = 1.0e-6;
	    }
	    else if (e > (1.0 - 1.0e-6))
	    {
	        e = 1.0 - 1.0e-6;
	    }

	    /*
	     * using calculated values, find position and velocity
	     * we can pass in finalants from Initialise() as these dont change
	     */
	    return CalculateFinalPositionVelocity(tsince, e,
	            a, omega, xl, xnode,
	            xincl, common_consts_.xlcof, common_consts_.aycof,
	            common_consts_.x3thm1, common_consts_.x1mth2, common_consts_.x7thm1,
	            common_consts_.cosio, common_consts_.sinio);

	}

	Eci CalculateFinalPositionVelocity(
	        final double tsince,
	        final double e,
	        final double a,
	        final double omega,
	        final double xl,
	        final double xnode,
	        final double xincl,
	        final double xlcof,
	        final double aycof,
	        final double x3thm1,
	        final double x1mth2,
	        final double x7thm1,
	        final double cosio,
	        final double sinio) throws SatelliteException, DecayedException
	{
	    final double beta2 = 1.0 - e * e;
	    final double xn = Globals.kXKE / Math.pow(a, 1.5);
	    /*
	     * long period periodics
	     */
	    final double axn = e * Math.cos(omega);
	    final double temp11 = 1.0 / (a * beta2);
	    final double xll = temp11 * xlcof * axn;
	    final double aynl = temp11 * aycof;
	    final double xlt = xl + xll;
	    final double ayn = e * Math.sin(omega) + aynl;
	    final double elsq = axn * axn + ayn * ayn;

	    if (elsq >= 1.0)
	    {
	        throw new SatelliteException("Error: (elsq >= 1.0)");
	    }

	    /*
	     * solve keplers equation
	     * - solve using Newton-Raphson root solving
	     * - here capu is almost the mean anomoly
	     * - initialise the eccentric anomaly term epw
	     * - The fmod saves reduction of angle to +/-2pi in sin/cos() and prevents
	     * convergence problems.
	     */
	    final double capu = (xlt - xnode) % Globals.kTWOPI;
	    double epw = capu;

	    double sinepw = 0.0;
	    double cosepw = 0.0;
	    double ecose = 0.0;
	    double esine = 0.0;

	    /*
	     * sensibility check for N-R correction
	     */
	    final double max_newton_naphson = 1.25 * Math.abs(Math.sqrt(elsq));

	    boolean kepler_running = true;

	    for (int i = 0; i < 10 && kepler_running; i++)
	    {
	        sinepw = Math.sin(epw);
	        cosepw = Math.cos(epw);
	        ecose = axn * cosepw + ayn * sinepw;
	        esine = axn * sinepw - ayn * cosepw;

	        double f = capu - epw + esine;

	        if (Math.abs(f) < 1.0e-12)
	        {
	            kepler_running = false;
	        }
	        else
	        {
	            /*
	             * 1st order Newton-Raphson correction
	             */
	            final double fdot = 1.0 - ecose;
	            double delta_epw = f / fdot;

	            /*
	             * 2nd order Newton-Raphson correction.
	             * f / (fdot - 0.5 * d2f * f/fdot)
	             */
	            if (i == 0)
	            {
	                if (delta_epw > max_newton_naphson)
	                {
	                    delta_epw = max_newton_naphson;
	                }
	                else if (delta_epw < -max_newton_naphson)
	                {
	                    delta_epw = -max_newton_naphson;
	                }
	            }
	            else
	            {
	                delta_epw = f / (fdot + 0.5 * esine * delta_epw);
	            }

	            /*
	             * Newton-Raphson correction of -F/DF
	             */
	            epw += delta_epw;
	        }
	    }
	    /*
	     * short period preliminary quantities
	     */
	    final double temp21 = 1.0 - elsq;
	    final double pl = a * temp21;

	    if (pl < 0.0)
	    {
	        throw new SatelliteException("Error: (pl < 0.0)");
	    }

	    final double r = a * (1.0 - ecose);
	    final double temp31 = 1.0 / r;
	    final double rdot = Globals.kXKE * Math.sqrt(a) * esine * temp31;
	    final double rfdot = Globals.kXKE * Math.sqrt(pl) * temp31;
	    final double temp32 = a * temp31;
	    final double betal = Math.sqrt(temp21);
	    final double temp33 = 1.0 / (1.0 + betal);
	    final double cosu = temp32 * (cosepw - axn + ayn * esine * temp33);
	    final double sinu = temp32 * (sinepw - ayn - axn * esine * temp33);
	    final double u = Math.atan2(sinu, cosu);
	    final double sin2u = 2.0 * sinu * cosu;
	    final double cos2u = 2.0 * cosu * cosu - 1.0;

	    /*
	     * update for short periodics
	     */
	    final double temp41 = 1.0 / pl;
	    final double temp42 = Globals.kCK2 * temp41;
	    final double temp43 = temp42 * temp41;

	    final double rk = r * (1.0 - 1.5 * temp43 * betal * x3thm1)
	        + 0.5 * temp42 * x1mth2 * cos2u;
	    final double uk = u - 0.25 * temp43 * x7thm1 * sin2u;
	    final double xnodek = xnode + 1.5 * temp43 * cosio * sin2u;
	    final double xinck = xincl + 1.5 * temp43 * cosio * sinio * cos2u;
	    final double rdotk = rdot - xn * temp42 * x1mth2 * sin2u;
	    final double rfdotk = rfdot + xn * temp42 * (x1mth2 * cos2u + 1.5 * x3thm1);

	    /*
	     * orientation vectors
	     */
	    final double sinuk = Math.sin(uk);
	    final double cosuk = Math.cos(uk);
	    final double sinik = Math.sin(xinck);
	    final double cosik = Math.cos(xinck);
	    final double sinnok = Math.sin(xnodek);
	    final double cosnok = Math.cos(xnodek);
	    final double xmx = -sinnok * cosik;
	    final double xmy = cosnok * cosik;
	    final double ux = xmx * sinuk + cosnok * cosuk;
	    final double uy = xmy * sinuk + sinnok * cosuk;
	    final double uz = sinik * sinuk;
	    final double vx = xmx * cosuk - cosnok * sinuk;
	    final double vy = xmy * cosuk - sinnok * sinuk;
	    final double vz = sinik * cosuk;
	    /*
	     * position and velocity
	     */
	    final double x = rk * ux * Globals.kXKMPER;
	    final double y = rk * uy * Globals.kXKMPER;
	    final double z = rk * uz * Globals.kXKMPER;
	    Vector4d position = new Vector4d(x, y, z,0);
	    final double xdot = (rdotk * ux + rfdotk * vx) * Globals.kXKMPER / 60.0;
	    final double ydot = (rdotk * uy + rfdotk * vy) * Globals.kXKMPER / 60.0;
	    final double zdot = (rdotk * uz + rfdotk * vz) * Globals.kXKMPER / 60.0;
	    Vector4d velocity= new Vector4d(xdot, ydot, zdot,0);

	    LocalDateTime dd = elements_.Epoch().plus((long)(tsince*1000), ChronoUnit.MILLIS);
	    if (rk < 1.0)
	    {
	        throw new DecayedException( dd, position, velocity);
	    }
	    
	    
	    return new Eci(dd, position, velocity);
	}

	static double EvaluateCubicPolynomial(
	        final double x,
	        final double finalant,
	        final double linear,
	        final double squared,
	        final double cubed)
	{
	    return finalant + x * (linear + x * (squared + x * cubed));
	}

	void DeepSpaceInitialise(
	        final double eosq,
	        final double sinio,
	        final double cosio,
	        final double betao,
	        final double theta2,
	        final double betao2,
	        final double xmdot,
	        final double omgdot,
	        final double xnodot)
	{
	    double se = 0.0;
	    double si = 0.0;
	    double sl = 0.0;
	    double sgh = 0.0;
	    double shdq = 0.0;

	    double bfact = 0.0;



	    final double aqnv = 1.0 / elements_.RecoveredSemiMajorAxis();
	    final double xpidot = omgdot + xnodot;
	    final double sinq = Math.sin(elements_.AscendingNode());
	    final double cosq = Math.cos(elements_.AscendingNode());
	    final double sing = Math.sin(elements_.ArgumentPerigee());
	    final double cosg = Math.cos(elements_.ArgumentPerigee());

	    /*
	     * initialize lunar / solar terms
	     */
	    double d =Utils.toJulianDay(Utils.getDoubleUnixTime(elements_.Epoch()));//todo: check
	    final double jday = (d - Globals.kEPOCH_JAN1_12H_2000);
	    
	    final double xnodce = 4.5236020 - 9.2422029e-4 * jday;
	    final double xnodce_temp = (xnodce% Globals.kTWOPI);
	    final double stem = Math.sin(xnodce_temp);
	    final double ctem = Math.cos(xnodce_temp);
	    final double zcosil = 0.91375164 - 0.03568096 * ctem;
	    final double zsinil = Math.sqrt(1.0 - zcosil * zcosil);
	    final double zsinhl = 0.089683511 * stem / zsinil;
	    final double zcoshl = Math.sqrt(1.0 - zsinhl * zsinhl);
	    final double c = 4.7199672 + 0.22997150 * jday;
	    final double gam = 5.8351514 + 0.0019443680 * jday;
	    deepspace_consts_.zmol = Utils.WrapTwoPI(c - gam);
	    double zx = 0.39785416 * stem / zsinil;
	    double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
	    zx = Math.atan2(zx, zy);
	    zx = (gam + zx - xnodce) % Globals.kTWOPI;

	    final double zcosgl = Math.cos(zx);
	    final double zsingl = Math.sin(zx);
	    deepspace_consts_.zmos = Utils.WrapTwoPI(6.2565837 + 0.017201977 * jday);

	    /*
	     * do solar terms
	     */
	    double zcosg = ZCOSGS;
	    double zsing = ZSINGS;
	    double zcosi = ZCOSIS;
	    double zsini = ZSINI;
	    double zcosh = cosq;
	    double zsinh = sinq;
	    double cc = C1SS;
	    double zn = ZNS;
	    double ze = ZES;
	    final double xnoi = 1.0 / elements_.RecoveredMeanMotion();

	    for (int cnt = 0; cnt < 2; cnt++)
	    {
	        /*
	         * solar terms are done a second time after lunar terms are done
	         */
	        final double a1 = zcosg * zcosh + zsing * zcosi * zsinh;
	        final double a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
	        final double a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
	        final double a8 = zsing * zsini;
	        final double a9 = zsing * zsinh + zcosg * zcosi*zcosh;
	        final double a10 = zcosg * zsini;
	        final double a2 = cosio * a7 + sinio * a8;
	        final double a4 = cosio * a9 + sinio * a10;
	        final double a5 = -sinio * a7 + cosio * a8;
	        final double a6 = -sinio * a9 + cosio * a10;
	        final double x1 = a1 * cosg + a2 * sing;
	        final double x2 = a3 * cosg + a4 * sing;
	        final double x3 = -a1 * sing + a2 * cosg;
	        final double x4 = -a3 * sing + a4 * cosg;
	        final double x5 = a5 * sing;
	        final double x6 = a6 * sing;
	        final double x7 = a5 * cosg;
	        final double x8 = a6 * cosg;
	        final double z31 = 12.0 * x1 * x1 - 3. * x3 * x3;
	        final double z32 = 24.0 * x1 * x2 - 6. * x3 * x4;
	        final double z33 = 12.0 * x2 * x2 - 3. * x4 * x4;
	        double z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * eosq;
	        double z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * eosq;
	        double z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * eosq;

	        final double z11 = -6.0 * a1 * a5
	            + eosq * (-24. * x1 * x7 - 6. * x3 * x5);
	        final double z12 = -6.0 * (a1 * a6 + a3 * a5) 
	            + eosq * (-24. * (x2 * x7 + x1 * x8) - 6. * (x3 * x6 + x4 * x5));
	        final double z13 = -6.0 * a3 * a6
	            + eosq * (-24. * x2 * x8 - 6. * x4 * x6);
	        final double z21 = 6.0 * a2 * a5
	            + eosq * (24. * x1 * x5 - 6. * x3 * x7);
	        final double z22 = 6.0 * (a4 * a5 + a2 * a6)
	            + eosq * (24. * (x2 * x5 + x1 * x6) - 6. * (x4 * x7 + x3 * x8));
	        final double z23 = 6.0 * a4 * a6
	            + eosq * (24. * x2 * x6 - 6. * x4 * x8);

	        z1 = z1 + z1 + betao2 * z31;
	        z2 = z2 + z2 + betao2 * z32;
	        z3 = z3 + z3 + betao2 * z33;

	        final double s3 = cc * xnoi;
	        final double s2 = -0.5 * s3 / betao;
	        final double s4 = s3 * betao;
	        final double s1 = -15.0 * elements_.Eccentricity() * s4;
	        final double s5 = x1 * x3 + x2 * x4;
	        final double s6 = x2 * x3 + x1 * x4;
	        final double s7 = x2 * x4 - x1 * x3;

	        se = s1 * zn * s5;
	        si = s2 * zn * (z11 + z13);
	        sl = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eosq);
	        sgh = s4 * zn * (z31 + z33 - 6.0);

	        /*
	         * replaced
	         * sh = -zn * s2 * (z21 + z23
	         * with
	         * shdq = (-zn * s2 * (z21 + z23)) / sinio
	         */
	        if (elements_.Inclination() < 5.2359877e-2
	                || elements_.Inclination() > Globals.kPI - 5.2359877e-2)
	        {
	            shdq = 0.0;
	        }
	        else
	        {
	            shdq = (-zn * s2 * (z21 + z23)) / sinio;
	        }

	        deepspace_consts_.ee2 = 2.0 * s1 * s6;
	        deepspace_consts_.e3 = 2.0 * s1 * s7;
	        deepspace_consts_.xi2 = 2.0 * s2 * z12;
	        deepspace_consts_.xi3 = 2.0 * s2 * (z13 - z11);
	        deepspace_consts_.xl2 = -2.0 * s3 * z2;
	        deepspace_consts_.xl3 = -2.0 * s3 * (z3 - z1);
	        deepspace_consts_.xl4 = -2.0 * s3 * (-21.0 - 9.0 * eosq) * ze;
	        deepspace_consts_.xgh2 = 2.0 * s4 * z32;
	        deepspace_consts_.xgh3 = 2.0 * s4 * (z33 - z31);
	        deepspace_consts_.xgh4 = -18.0 * s4 * ze;
	        deepspace_consts_.xh2 = -2.0 * s2 * z22;
	        deepspace_consts_.xh3 = -2.0 * s2 * (z23 - z21);

	        if (cnt == 1)
	        {
	            break;
	        }
	        /*
	         * do lunar terms
	         */
	        deepspace_consts_.sse = se;
	        deepspace_consts_.ssi = si;
	        deepspace_consts_.ssl = sl;
	        deepspace_consts_.ssh = shdq;
	        deepspace_consts_.ssg = sgh - cosio * deepspace_consts_.ssh;
	        deepspace_consts_.se2 = deepspace_consts_.ee2;
	        deepspace_consts_.si2 = deepspace_consts_.xi2;
	        deepspace_consts_.sl2 = deepspace_consts_.xl2;
	        deepspace_consts_.sgh2 = deepspace_consts_.xgh2;
	        deepspace_consts_.sh2 = deepspace_consts_.xh2;
	        deepspace_consts_.se3 = deepspace_consts_.e3;
	        deepspace_consts_.si3 = deepspace_consts_.xi3;
	        deepspace_consts_.sl3 = deepspace_consts_.xl3;
	        deepspace_consts_.sgh3 = deepspace_consts_.xgh3;
	        deepspace_consts_.sh3 = deepspace_consts_.xh3;
	        deepspace_consts_.sl4 = deepspace_consts_.xl4;
	        deepspace_consts_.sgh4 = deepspace_consts_.xgh4;
	        zcosg = zcosgl;
	        zsing = zsingl;
	        zcosi = zcosil;
	        zsini = zsinil;
	        zcosh = zcoshl * cosq + zsinhl * sinq;
	        zsinh = sinq * zcoshl - cosq * zsinhl;
	        zn = ZNL;
	        cc = C1L;
	        ze = ZEL;
	    }

	    deepspace_consts_.sse += se;
	    deepspace_consts_.ssi += si;
	    deepspace_consts_.ssl += sl;
	    deepspace_consts_.ssg += sgh - cosio * shdq;
	    deepspace_consts_.ssh += shdq;

	    deepspace_consts_.resonance_flag = false;
	    deepspace_consts_.synchronous_flag = false;
	    boolean initialise_integrator = true;

	    if (elements_.RecoveredMeanMotion() < 0.0052359877
	            && elements_.RecoveredMeanMotion() > 0.0034906585)
	    {
	        /*
	         * 24h synchronous resonance terms initialisation
	         */
	        deepspace_consts_.resonance_flag = true;
	        deepspace_consts_.synchronous_flag = true;

	        final double g200 = 1.0 + eosq * (-2.5 + 0.8125 * eosq);
	        final double g310 = 1.0 + 2.0 * eosq;
	        final double g300 = 1.0 + eosq * (-6.0 + 6.60937 * eosq);
	        final double f220 = 0.75 * (1.0 + cosio) * (1.0 + cosio);
	        final double f311 = 0.9375 * sinio * sinio * (1.0 + 3.0 * cosio)
	            - 0.75 * (1.0 + cosio);
	        double f330 = 1.0 + cosio;
	        f330 = 1.875 * f330 * f330 * f330;
	        deepspace_consts_.del1 = 3.0 * elements_.RecoveredMeanMotion()
	            * elements_.RecoveredMeanMotion()
	            * aqnv * aqnv;
	        deepspace_consts_.del2 = 2.0 * deepspace_consts_.del1
	            * f220 * g200 * Q22;
	        deepspace_consts_.del3 = 3.0 * deepspace_consts_.del1
	            * f330 * g300 * Q33 * aqnv;
	        deepspace_consts_.del1 = deepspace_consts_.del1
	            * f311 * g310 * Q31 * aqnv;

	        integrator_consts_.xlamo = elements_.MeanAnomoly()
	            + elements_.AscendingNode()
	            + elements_.ArgumentPerigee()
	            - deepspace_consts_.gsto;
	        bfact = xmdot + xpidot - Globals.kTHDT;
	        bfact += deepspace_consts_.ssl
	            + deepspace_consts_.ssg
	            + deepspace_consts_.ssh;
	    }
	    else if (elements_.RecoveredMeanMotion() < 8.26e-3
	            || elements_.RecoveredMeanMotion() > 9.24e-3
	            || elements_.Eccentricity() < 0.5)
	    {
	        initialise_integrator = false;
	    }
	    else
	    {
	        /*
	         * geopotential resonance initialisation for 12 hour orbits
	         */
	        deepspace_consts_.resonance_flag = true;

	        double g211;
	        double g310;
	        double g322;
	        double g410;
	        double g422;
	        double g520;

	        double g201 = -0.306 - (elements_.Eccentricity() - 0.64) * 0.440;

	        if (elements_.Eccentricity() <= 0.65)
	        {
	            g211 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    3.616, -13.247, +16.290, 0.0);
	            g310 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -19.302, 117.390, -228.419, 156.591);
	            g322 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -18.9068, 109.7927, -214.6334, 146.5816);
	            g410 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -41.122, 242.694, -471.094, 313.953);
	            g422 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -146.407, 841.880, -1629.014, 1083.435);
	            g520 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -532.114, 3017.977, -5740.032, 3708.276);
	        }
	        else
	        {
	            g211 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -72.099, 331.819, -508.738, 266.724);
	            g310 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -346.844, 1582.851, -2415.925, 1246.113);
	            g322 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -342.585, 1554.908, -2366.899, 1215.972);
	            g410 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -1052.797, 4758.686, -7193.992, 3651.957);
	            g422 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -3581.69, 16178.11, -24462.77, 12422.52);

	            if (elements_.Eccentricity() <= 0.715)
	            {
	                g520 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                        1464.74, -4664.75, 3763.64, 0.0);
	            }
	            else
	            {
	                g520 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                        -5149.66, 29936.92, -54087.36, 31324.56);
	            }
	        }

	        double g533;
	        double g521;
	        double g532;

	        if (elements_.Eccentricity() < 0.7)
	        {
	            g533 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -919.2277, 4988.61, -9064.77, 5542.21);
	            g521 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -822.71072, 4568.6173, -8491.4146, 5337.524);
	            g532 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -853.666, 4690.25, -8624.77, 5341.4);
	        }
	        else
	        {
	            g533 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -37995.78, 161616.52, -229838.2, 109377.94);
	            g521 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -51752.104, 218913.95, -309468.16, 146349.42);
	            g532 = EvaluateCubicPolynomial(elements_.Eccentricity(),
	                    -40023.88, 170470.89, -242699.48, 115605.82);
	        }

	        final double sini2 = sinio * sinio;
	        final double f220 = 0.75 * (1.0 + 2.0 * cosio + theta2);
	        final double f221 = 1.5 * sini2;
	        final double f321 = 1.875 * sinio * (1.0 - 2.0 * cosio - 3.0 * theta2);
	        final double f322 = -1.875 * sinio * (1.0 + 2.0 * cosio - 3.0 * theta2);
	        final double f441 = 35.0 * sini2 * f220;
	        final double f442 = 39.3750 * sini2 * sini2;
	        final double f522 = 9.84375 * sinio
	            * (sini2 * (1.0 - 2.0 * cosio - 5.0 * theta2)
	                + 0.33333333 * (-2.0 + 4.0 * cosio + 6.0 * theta2));
	        final double f523 = sinio
	            * (4.92187512 * sini2 * (-2.0 - 4.0 * cosio + 10.0 * theta2)
	                + 6.56250012 * (1.0 + 2.0 * cosio - 3.0 * theta2));
	        final double f542 = 29.53125 * sinio * (2.0 - 8.0 * cosio + theta2 *
	                (-12.0 + 8.0 * cosio + 10.0 * theta2));
	        final double f543 = 29.53125 * sinio * (-2.0 - 8.0 * cosio + theta2 *
	                (12.0 + 8.0 * cosio - 10.0 * theta2));

	        final double xno2 = elements_.RecoveredMeanMotion()
	            * elements_.RecoveredMeanMotion();
	        final double ainv2 = aqnv * aqnv;

	        double temp1 = 3.0 * xno2 * ainv2;
	        double temp = temp1 * ROOT22;
	        deepspace_consts_.d2201 = temp * f220 * g201;
	        deepspace_consts_.d2211 = temp * f221 * g211;
	        temp1 = temp1 * aqnv;
	        temp = temp1 * ROOT32;
	        deepspace_consts_.d3210 = temp * f321 * g310;
	        deepspace_consts_.d3222 = temp * f322 * g322;
	        temp1 = temp1 * aqnv;
	        temp = 2.0 * temp1 * ROOT44;
	        deepspace_consts_.d4410 = temp * f441 * g410;
	        deepspace_consts_.d4422 = temp * f442 * g422;
	        temp1 = temp1 * aqnv;
	        temp = temp1 * ROOT52;
	        deepspace_consts_.d5220 = temp * f522 * g520;
	        deepspace_consts_.d5232 = temp * f523 * g532;
	        temp = 2.0 * temp1 * ROOT54;
	        deepspace_consts_.d5421 = temp * f542 * g521;
	        deepspace_consts_.d5433 = temp * f543 * g533;

	        integrator_consts_.xlamo = elements_.MeanAnomoly()
	            + elements_.AscendingNode()
	            + elements_.AscendingNode()
	            - deepspace_consts_.gsto
	            - deepspace_consts_.gsto;
	        bfact = xmdot
	            + xnodot + xnodot
	            - Globals.kTHDT - Globals.kTHDT;
	        bfact = bfact + deepspace_consts_.ssl
	            + deepspace_consts_.ssh
	            + deepspace_consts_.ssh;
	    }

	    if (initialise_integrator)
	    {
	        /*
	         * initialise integrator
	         */
	        integrator_consts_.xfact = bfact - elements_.RecoveredMeanMotion();
	        integrator_params_.atime = 0.0;
	        integrator_params_.xni = elements_.RecoveredMeanMotion();
	        integrator_params_.xli = integrator_consts_.xlamo;
	        /*
	         * precompute dot terms for epoch
	         */
	        DeepSpaceCalcDotTerms(integrator_consts_.values_0);
	    }
	}

	void DeepSpaceCalculateLunarSolarTerms(final double tsince,double[] pe, double[] pinc, double[] pl,double[] pgh, double[] ph)
	{
	    //static final double ZES = 0.01675;
	    //static final double ZNS = 1.19459E-5;
	    //static final double ZNL = 1.5835218E-4;
	    //static final double ZEL = 0.05490;

	    /*
	     * calculate solar terms for time tsince
	     */
	    double zm = deepspace_consts_.zmos + ZNS * tsince;
	    double zf = zm + 2.0 * ZES * Math.sin(zm);
	    double sinzf = Math.sin(zf);
	    double f2 = 0.5 * sinzf * sinzf - 0.25;
	    double f3 = -0.5 * sinzf * Math.cos(zf);

	    final double ses = deepspace_consts_.se2 * f2
	        + deepspace_consts_.se3 * f3;
	    final double sis = deepspace_consts_.si2 * f2
	        + deepspace_consts_.si3 * f3;
	    final double sls = deepspace_consts_.sl2 * f2
	        + deepspace_consts_.sl3 * f3
	        + deepspace_consts_.sl4 * sinzf;
	    final double sghs = deepspace_consts_.sgh2 * f2
	        + deepspace_consts_.sgh3 * f3
	        + deepspace_consts_.sgh4 * sinzf;
	    final double shs = deepspace_consts_.sh2 * f2
	        + deepspace_consts_.sh3 * f3;

	    /*
	     * calculate lunar terms for time tsince
	     */
	    zm = deepspace_consts_.zmol + ZNL * tsince;
	    zf = zm + 2.0 * ZEL * Math.sin(zm);
	    sinzf = Math.sin(zf);
	    f2 = 0.5 * sinzf * sinzf - 0.25;
	    f3 = -0.5 * sinzf * Math.cos(zf);

	    final double sel = deepspace_consts_.ee2 * f2
	        + deepspace_consts_.e3 * f3;
	    final double sil = deepspace_consts_.xi2 * f2
	        + deepspace_consts_.xi3 * f3;
	    final double sll = deepspace_consts_.xl2 * f2
	        + deepspace_consts_.xl3 * f3
	        + deepspace_consts_.xl4 * sinzf;
	    final double sghl = deepspace_consts_.xgh2 * f2
	        + deepspace_consts_.xgh3 * f3
	        + deepspace_consts_.xgh4 * sinzf;
	    final double shl = deepspace_consts_.xh2 * f2
	        + deepspace_consts_.xh3 * f3;

	    /*
	     * merge calculated values
	     */
	    pe[0] = ses + sel;
	    pinc[0] = sis + sil;
	    pl[0] = sls + sll;
	    pgh[0] = sghs + sghl;
	    ph[0] = shs + shl;
	}

	void DeepSpacePeriodics(
	        final double tsince,
	        double[] em,
	        double[] xinc,
	        double[] omgasm,
	        double[] xnodes,
	        double[] xll)
	{
	    /*
	     * storage for lunar / solar terms
	     * set by DeepSpaceCalculateLunarSolarTerms()
	     */
	    double[] pe = new double[]{0.0};
	    double[] pinc = new double[]{0.0};
	    double[] pl = new double[]{0.0};
	    double[] pgh = new double[]{0.0};
	    double[] ph = new double[]{0.0};

	    /*
	     * calculate lunar / solar terms for current time
	     */
	    DeepSpaceCalculateLunarSolarTerms(tsince, pe, pinc, pl, pgh, ph);

	    xinc[0] += pinc[0];
	    em[0] += pe[0];

	    /* Spacetrack report #3 has sin/cos from before perturbations
	     * added to xinc (oldxinc), but apparently report # 6 has then
	     * from after they are added.
	     * use for strn3
	     * if (elements_.Inclination() >= 0.2)
	     * use for gsfc
	     * if (xinc >= 0.2)
	     * (moved from start of function)
	     */
	    final double sinis = Math.sin(xinc[0]);
	    final double cosis = Math.cos(xinc[0]);

	    if (xinc[0] >= 0.2)
	    {
	        /*
	         * apply periodics directly
	         */
	        final double tmp_ph = ph[0] / sinis;

	        omgasm[0] += pgh[0] - cosis * tmp_ph;
	        xnodes[0] += tmp_ph;
	        xll[0] += pl[0];
	    }
	    else
	    {
	        /*
	         * apply periodics with lyddane modification
	         */
	        final double sinok = Math.sin(xnodes[0]);
	        final double cosok = Math.cos(xnodes[0]);
	        double alfdp = sinis * sinok;
	        double betdp = sinis * cosok;
	        final double dalf = ph[0] * cosok + pinc[0] * cosis * sinok;
	        final double dbet = -ph[0] * sinok + pinc[0] * cosis * cosok;

	        alfdp += dalf;
	        betdp += dbet;

	        xnodes[0] = Utils.WrapTwoPI(xnodes[0]);

	        double xls = xll[0] + omgasm[0] + cosis * xnodes[0];
	        double dls = pl[0] + pgh[0] - pinc[0] * xnodes[0] * sinis;
	        xls += dls;

	        /*
	         * save old xnodes value
	         */
	        final double oldxnodes = xnodes[0];

	        xnodes[0] = Math.atan2(alfdp, betdp);
	        if (xnodes[0] < 0.0)
	        {
	            xnodes[0] += Globals.kTWOPI;
	        }

	        /*
	         * Get perturbed xnodes in to same quadrant as original.
	         * RAAN is in the range of 0 to 360 degrees
	         * atan2 is in the range of -180 to 180 degrees
	         */
	        if (Math.abs(oldxnodes - xnodes[0]) > Globals.kPI)
	        {
	            if (xnodes[0] < oldxnodes)
	            {
	                xnodes[0] += Globals.kTWOPI;
	            }
	            else
	            {
	                xnodes[0] -= Globals.kTWOPI;
	            }
	        }

	        xll[0] += pl[0];
	        omgasm[0] = xls - xll[0] - cosis * xnodes[0];
	    }
	}
	    static final double STEP = 720.0;
	    static final double STEP2 = 259200.0;
	void DeepSpaceSecular(
	        final double tsince,
	        double[] xll,
	        double[] omgasm,
	        double[] xnodes,
	        double[] em,
	        double[] xinc,
	        double[] xn) 
	{


	    xll[0] += deepspace_consts_.ssl * tsince;
	    omgasm[0] += deepspace_consts_.ssg * tsince;
	    xnodes[0] += deepspace_consts_.ssh * tsince;
	    em[0] += deepspace_consts_.sse * tsince;
	    xinc[0] += deepspace_consts_.ssi * tsince;

	    if (deepspace_consts_.resonance_flag)
	    {
	        /*
	         * 1st condition (if tsince is less than one time step from epoch)
	         * 2nd condition (if integrator_params_.atime and
	         *     tsince are of opposite signs, so zero crossing required)
	         * 3rd condition (if tsince is closer to zero than 
	         *     integrator_params_.atime, only integrate away from zero)
	         */
	        if (Math.abs(tsince) < STEP ||
	                tsince * integrator_params_.atime <= 0.0 ||
	                Math.abs(tsince) < Math.abs(integrator_params_.atime))
	        {
	            /*
	             * restart from epoch
	             */
	            integrator_params_.atime = 0.0;
	            integrator_params_.xni = elements_.RecoveredMeanMotion();
	            integrator_params_.xli = integrator_consts_.xlamo;

	            /*
	             * restore precomputed values for epoch
	             */
	            integrator_params_.values_t = integrator_consts_.values_0;
	        }

	        double ft = tsince - integrator_params_.atime;

	        /*
	         * if time difference (ft) is greater than the time step (720.0)
	         * loop around until integrator_params_.atime is within one time step of
	         * tsince
	         */
	        if (Math.abs(ft) >= STEP)
	        {
	            /*
	             * calculate step direction to allow integrator_params_.atime
	             * to catch up with tsince
	             */
	            double delt = -STEP;
	            if (ft >= 0.0)
	            {
	                delt = STEP;
	            }

	            do
	            {
	                /*
	                 * integrate using current dot terms
	                 */
	                DeepSpaceIntegrator(delt, STEP2, integrator_params_.values_t);

	                /*
	                 * calculate dot terms for next integration
	                 */
	                DeepSpaceCalcDotTerms(integrator_params_.values_t);

	                ft = tsince - integrator_params_.atime;
	            } while (Math.abs(ft) >= STEP);
	        }

	        /*
	         * integrator
	         */
	        xn[0] = integrator_params_.xni 
	            + integrator_params_.values_t.xndot * ft
	            + integrator_params_.values_t.xnddt * ft * ft * 0.5;
	        final double xl = integrator_params_.xli
	            + integrator_params_.values_t.xldot * ft
	            + integrator_params_.values_t.xndot * ft * ft * 0.5;
	        final double temp = -xnodes[0] + deepspace_consts_.gsto + tsince * Globals.kTHDT;

	        if (deepspace_consts_.synchronous_flag)
	        {
	            xll[0] = xl + temp - omgasm[0];
	        }
	        else
	        {
	            xll[0] = xl + temp + temp;
	        }
	    }
	}
    static final double G22 = 5.7686396;
    static final double G32 = 0.95240898;
    static final double G44 = 1.8014998;
    static final double G52 = 1.0508330;
    static final double G54 = 4.4108898;
    static final double FASX2 = 0.13130908;
    static final double FASX4 = 2.8843198;
    static final double FASX6 = 0.37448087;
	void DeepSpaceCalcDotTerms(IntegratorValues values) 
	{


	    if (deepspace_consts_.synchronous_flag)
	    {

	        values.xndot = deepspace_consts_.del1
	            * Math.sin(integrator_params_.xli - FASX2)
	            + deepspace_consts_.del2
	            * Math.sin(2.0 * (integrator_params_.xli - FASX4))
	            + deepspace_consts_.del3
	            * Math.sin(3.0 * (integrator_params_.xli - FASX6));
	        values.xnddt = deepspace_consts_.del1
	            * Math.cos(integrator_params_.xli - FASX2)
	            + 2.0 * deepspace_consts_.del2
	            * Math.cos(2.0 * (integrator_params_.xli - FASX4))
	            + 3.0 * deepspace_consts_.del3
	            * Math.cos(3.0 * (integrator_params_.xli - FASX6));
	    }
	    else
	    {
	        final double xomi = elements_.ArgumentPerigee()
	            + common_consts_.omgdot * integrator_params_.atime;
	        final double x2omi = xomi + xomi;
	        final double x2li = integrator_params_.xli + integrator_params_.xli;

	        values.xndot = deepspace_consts_.d2201
	            * Math.sin(x2omi + integrator_params_.xli - G22)
	            * + deepspace_consts_.d2211
	            * Math.sin(integrator_params_.xli - G22)
	            + deepspace_consts_.d3210
	            * Math.sin(xomi + integrator_params_.xli - G32)
	            + deepspace_consts_.d3222
	            * Math.sin(-xomi + integrator_params_.xli - G32)
	            + deepspace_consts_.d4410
	            * Math.sin(x2omi + x2li - G44)
	            + deepspace_consts_.d4422
	            * Math.sin(x2li - G44)
	            + deepspace_consts_.d5220
	            * Math.sin(xomi + integrator_params_.xli - G52)
	            + deepspace_consts_.d5232
	            * Math.sin(-xomi + integrator_params_.xli - G52)
	            + deepspace_consts_.d5421
	            * Math.sin(xomi + x2li - G54)
	            + deepspace_consts_.d5433
	            * Math.sin(-xomi + x2li - G54);
	        values.xnddt = deepspace_consts_.d2201
	            * Math.cos(x2omi + integrator_params_.xli - G22)
	            + deepspace_consts_.d2211
	            * Math.cos(integrator_params_.xli - G22)
	            + deepspace_consts_.d3210
	            * Math.cos(xomi + integrator_params_.xli - G32)
	            + deepspace_consts_.d3222
	            * Math.cos(-xomi + integrator_params_.xli - G32)
	            + deepspace_consts_.d5220
	            * Math.cos(xomi + integrator_params_.xli - G52)
	            + deepspace_consts_.d5232
	            * Math.cos(-xomi + integrator_params_.xli - G52)
	            + 2.0 * (deepspace_consts_.d4410 * Math.cos(x2omi + x2li - G44)
	            + deepspace_consts_.d4422
	            * Math.cos(x2li - G44)
	            + deepspace_consts_.d5421
	            * Math.cos(xomi + x2li - G54)
	            + deepspace_consts_.d5433
	            * Math.cos(-xomi + x2li - G54));
	    }

	    values.xldot = integrator_params_.xni + integrator_consts_.xfact;
	    values.xnddt *= values.xldot;
	}

	void DeepSpaceIntegrator(
	        final double delt,
	        final double step2,
	        final IntegratorValues values)
	{
	    /*
	     * integrator
	     */
	    integrator_params_.xli += values.xldot * delt + values.xndot * step2;
	    integrator_params_.xni += values.xndot * delt + values.xnddt * step2;

	    /*
	     * increment integrator time
	     */
	    integrator_params_.atime += delt;
	}

	void Reset()
	{
	    use_simple_model_ = false;
	    use_deep_space_ = false;

	    common_consts_     = Empty_CommonConstants;
	    nearspace_consts_  = Empty_NearSpaceConstants;
	    deepspace_consts_  = Empty_DeepSpaceConstants;
	    integrator_consts_ = Empty_IntegratorConstants;
	    integrator_params_ = Empty_IntegratorParams;
	}

}
