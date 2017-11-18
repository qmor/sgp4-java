package org.vniiem.sgp4;
import java.time.LocalDateTime;



public class Tle
{

	/**
	 * @throws TleException 
	 * @details Initialise given the two lines of a tle
	 * @param[in] line_one Tle line one
	 * @param[in] line_two Tle line two
	 */
	public Tle(String line_one,   String line_two) throws TleException
	{
		line_one_ = line_one;
		line_two_ = line_two;
		Initialize();
	}

	/**
	 * @throws TleException 
	 * @details Initialise given the satellite name and the two lines of a tle
	 * @param[in] name Satellite name
	 * @param[in] line_one Tle line one
	 * @param[in] line_two Tle line two
	 */
	public Tle(String name,    		String line_one,    		String line_two) throws TleException
	{
		name_ = name;
		line_one_ = line_one;
		line_two_ = line_two;
		Initialize();
	}

	/**
	 * Copy constructor
	 * @param[in] tle Tle object to copy from
	 */
	Tle(Tle tle)
	{
		name_ = tle.name_;
		line_one_ = tle.line_one_;
		line_two_ = tle.line_two_;

		norad_number_ = tle.norad_number_;
		int_designator_ = tle.int_designator_;
		epoch_ = tle.epoch_;
		mean_motion_dt2_ = tle.mean_motion_dt2_;
		mean_motion_ddt6_ = tle.mean_motion_ddt6_;
		bstar_ = tle.bstar_;
		inclination_ = tle.inclination_;
		right_ascending_node_ = tle.right_ascending_node_;
		eccentricity_ = tle.eccentricity_;
		argument_perigee_ = tle.argument_perigee_;
		mean_anomaly_ = tle.mean_anomaly_;
		mean_motion_ = tle.mean_motion_;
		orbit_number_ = tle.orbit_number_;
	}

	/**
	 * Get the satellite name
	 * @returns the satellite name
	 */
	String Name() 
	{
		return name_;
	}

	/**
	 * Get the first line of the tle
	 * @returns the first line of the tle
	 */
	String Line1() 
	{
		return line_one_;
	}

	/**
	 * Get the second line of the tle
	 * @returns the second line of the tle
	 */
	String Line2() 
	{
		return line_two_;
	}

	/**
	 * Get the norad number
	 * @returns the norad number
	 */
	int NoradNumber()
	{
		return norad_number_;
	}

	/**
	 * Get the international designator
	 * @returns the international designator
	 */
	String IntDesignator()
	{
		return int_designator_;
	}

	/**
	 * Get the tle epoch
	 * @returns the tle epoch
	 */
	public LocalDateTime Epoch()
	{
		return epoch_;
	}

	/**
	 * Get the first time derivative of the mean motion divided by two
	 * @returns the first time derivative of the mean motion divided by two
	 */
	double MeanMotionDt2()
	{
		return mean_motion_dt2_;
	}

	/**
	 * Get the second time derivative of mean motion divided by six
	 * @returns the second time derivative of mean motion divided by six
	 */
	double MeanMotionDdt6()
	{
		return mean_motion_ddt6_;
	}

	/**
	 * Get the BSTAR drag term
	 * @returns the BSTAR drag term
	 */
	double BStar()
	{
		return bstar_;
	}

	/**
	 * Get the inclination
	 * @param in_degrees Whether to return the value in degrees or radians
	 * @returns the inclination
	 */
	double Inclination(Boolean in_degrees)
	{
		if (in_degrees)
		{
			return inclination_;
		}
		else
		{
			return Utils.DegreesToRadians(inclination_);
		}
	}

	/**
	 * Get the right ascension of the ascending node
	 * @param in_degrees Whether to return the value in degrees or radians
	 * @returns the right ascension of the ascending node
	 */
	double RightAscendingNode(Boolean in_degrees)
	{
		if (in_degrees)
		{
			return right_ascending_node_;
		}
		else
		{
			return Utils.DegreesToRadians(right_ascending_node_);
		}
	}

	/**
	 * Get the eccentricity
	 * @returns the eccentricity
	 */
	double Eccentricity() 
	{
		return eccentricity_;
	}

	/**
	 * Get the argument of perigee
	 * @param in_degrees Whether to return the value in degrees or radians
	 * @returns the argument of perigee
	 */
	double ArgumentPerigee(Boolean in_degrees) 
	{
		if (in_degrees)
		{
			return argument_perigee_;
		}
		else
		{
			return Utils.DegreesToRadians(argument_perigee_);
		}
	}

	/**
	 * Get the mean anomaly
	 * @param in_degrees Whether to return the value in degrees or radians
	 * @returns the mean anomaly
	 */
	double MeanAnomaly(Boolean in_degrees)
	{
		if (in_degrees)
		{
			return mean_anomaly_;
		}
		else
		{
			return Utils.DegreesToRadians(mean_anomaly_);
		}
	}

	/**
	 * Get the mean motion
	 * @returns the mean motion (revolutions per day)
	 */
	double MeanMotion()
	{
		return mean_motion_;
	}

	/**
	 * Get the orbit number
	 * @returns the orbit number
	 */
	int OrbitNumber()
	{
		return orbit_number_;
	}

	/**
	 * Get the expected tle line length
	 * @returns the tle line length
	 */
	static int LineLength()
	{
		return TLE_LEN_LINE_DATA;
	}

	/**
	 * Dump this object to a string
	 * @returns string
	 */
	@Override
	public
	String toString()
	{
		String ss="";


		ss+= "Norad Number:         " + NoradNumber()+"\n";
		ss += "Int. Designator:      " + IntDesignator() +"\n";
		ss += "Epoch:                " + Epoch() + "\n";
		ss  +="Orbit Number:         " + OrbitNumber() + "\n";
		ss  += "Mean Motion Dt2:      ";
		ss  += MeanMotionDt2() +"\n";
		ss  += "Mean Motion Ddt6:     ";
		ss  +=  MeanMotionDdt6() + "\n";
		ss  += "Eccentricity:         ";
		ss  +=  Eccentricity() + "\n";
		ss  += "BStar:                ";
		ss  += BStar() +"\n";
		ss  += "Inclination:          ";
		ss  += Inclination(true) + "\n";
		ss  += "Right Ascending Node: ";
		ss  += RightAscendingNode(true) + "\n";
		ss  += "Argument Perigee:     ";
		ss  += ArgumentPerigee(true) + "\n";
		ss  += "Mean Anomaly:         ";
		ss  +=  MeanAnomaly(true) + "\n";
		ss  +="Mean Motion:          ";
		ss  += MeanMotion() + "\n";
		return ss;
	}




	String name_;
	String line_one_;
	String line_two_;

	String int_designator_;
	LocalDateTime epoch_;
	double mean_motion_dt2_;
	double mean_motion_ddt6_;
	double bstar_;
	double inclination_;
	double right_ascending_node_;
	double eccentricity_;
	double argument_perigee_;
	double mean_anomaly_;
	double mean_motion_;
	int norad_number_;
	int orbit_number_;

	static final int TLE_LEN_LINE_DATA = 69;
	static final int TLE_LEN_LINE_NAME = 22;



	static final int TLE1_COL_NORADNUM = 2;
	static final int  TLE1_LEN_NORADNUM = 5;
	static final int TLE1_COL_INTLDESC_A = 9;
	static final int TLE1_LEN_INTLDESC_A = 2;
	//  static final int TLE1_COL_INTLDESC_B = 11;
	static final int TLE1_LEN_INTLDESC_B = 3;
	//  static final int TLE1_COL_INTLDESC_C = 14;
	static final int TLE1_LEN_INTLDESC_C = 3;
	static final int TLE1_COL_EPOCH_A = 18;
	static final int TLE1_LEN_EPOCH_A = 2;
	static final int TLE1_COL_EPOCH_B = 20;
	static final int TLE1_LEN_EPOCH_B = 12;
	static final int TLE1_COL_MEANMOTIONDT2 = 33;
	static final int TLE1_LEN_MEANMOTIONDT2 = 10;
	static final int TLE1_COL_MEANMOTIONDDT6 = 44;
	static final int TLE1_LEN_MEANMOTIONDDT6 = 8;
	static final int TLE1_COL_BSTAR = 53;
	static final int TLE1_LEN_BSTAR = 8;
	//  static final int TLE1_COL_EPHEMTYPE = 62;
	//  static final int TLE1_LEN_EPHEMTYPE = 1;
	//  static final int TLE1_COL_ELNUM = 64;
	//  static final int TLE1_LEN_ELNUM = 4;

	static final int TLE2_COL_NORADNUM = 2;
	static final int TLE2_LEN_NORADNUM = 5;
	static final int TLE2_COL_INCLINATION = 8;
	static final int TLE2_LEN_INCLINATION = 8;
	static final int TLE2_COL_RAASCENDNODE = 17;
	static final int TLE2_LEN_RAASCENDNODE = 8;
	static final int TLE2_COL_ECCENTRICITY = 26;
	static final int TLE2_LEN_ECCENTRICITY = 7;
	static final int TLE2_COL_ARGPERIGEE = 34;
	static final int TLE2_LEN_ARGPERIGEE = 8;
	static final int TLE2_COL_MEANANOMALY = 43;
	static final int TLE2_LEN_MEANANOMALY = 8;
	static final int TLE2_COL_MEANMOTION = 52;
	static final int TLE2_LEN_MEANMOTION = 11;
	static final int TLE2_COL_REVATEPOCH = 63;
	static final int TLE2_LEN_REVATEPOCH = 5;


	/**
	 * Initialise the tle object.
	 * @exception TleException
	 */
	void Initialize() throws TleException
	{
		if (!IsValidLineLength(line_one_))
		{
			throw new TleException("Invalid length for line one");
		}

		if (!IsValidLineLength(line_two_))
		{
			throw new TleException("Invalid length for line two");
		}

		if (line_one_.charAt(0) != '1')
		{
			throw new TleException("Invalid line beginning for line one");
		}

		if (line_two_.charAt(0) != '2')
		{
			throw new TleException("Invalid line beginning for line two");
		}

		int sat_number_1 = 0;
		int sat_number_2 = 0;

		sat_number_1 = ExtractInteger(Utils.SubString(line_one_,TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM));
		sat_number_2 =ExtractInteger(Utils.SubString(line_two_,TLE2_COL_NORADNUM, TLE2_LEN_NORADNUM));

		if (sat_number_1 != sat_number_2)
		{
			throw new TleException("Satellite numbers do not match");
		}

		norad_number_ = sat_number_1;

		if (name_.isEmpty())
		{
			name_ = line_one_.substring(TLE1_COL_NORADNUM, TLE1_LEN_NORADNUM);
		}

		int_designator_ = Utils.SubString(line_one_,TLE1_COL_INTLDESC_A,	TLE1_LEN_INTLDESC_A + TLE1_LEN_INTLDESC_B + TLE1_LEN_INTLDESC_C);

		int year = 0;
		double day = 0.0;

		year = ExtractInteger(Utils.SubString(line_one_,TLE1_COL_EPOCH_A,	TLE1_LEN_EPOCH_A));

		day = ExtractDouble(Utils.SubString(line_one_,TLE1_COL_EPOCH_B,TLE1_LEN_EPOCH_B), 4);
		mean_motion_dt2_ = ExtractDouble(Utils.SubString(line_one_,TLE1_COL_MEANMOTIONDT2,TLE1_LEN_MEANMOTIONDT2), 2);
		mean_motion_ddt6_= ExtractExponential(Utils.SubString(line_one_,TLE1_COL_MEANMOTIONDDT6, TLE1_LEN_MEANMOTIONDDT6));
		bstar_ = ExtractExponential(Utils.SubString(line_one_,TLE1_COL_BSTAR, TLE1_LEN_BSTAR));

		/*
		 * line 2
		 */
		inclination_ =  ExtractDouble(Utils.SubString(line_two_,TLE2_COL_INCLINATION, TLE2_LEN_INCLINATION), 4);
		right_ascending_node_ = ExtractDouble(Utils.SubString(line_two_,TLE2_COL_RAASCENDNODE, TLE2_LEN_RAASCENDNODE), 4);
		eccentricity_=  ExtractDouble(Utils.SubString(line_two_,TLE2_COL_ECCENTRICITY, TLE2_LEN_ECCENTRICITY), -1);
		argument_perigee_ = ExtractDouble(Utils.SubString(line_two_,TLE2_COL_ARGPERIGEE, TLE2_LEN_ARGPERIGEE), 4);
		mean_anomaly_ = ExtractDouble(Utils.SubString(line_two_,TLE2_COL_MEANANOMALY, TLE2_LEN_MEANANOMALY), 4);
		mean_motion_ = ExtractDouble(Utils.SubString(line_two_,TLE2_COL_MEANMOTION, TLE2_LEN_MEANMOTION), 3 );
		orbit_number_ =ExtractInteger(Utils.SubString(line_two_,TLE2_COL_REVATEPOCH, TLE2_LEN_REVATEPOCH));

		if (year < 57)
			year += 2000;
		else
			year += 1900;

		epoch_ =  Utils.epochDateTime(year, day);
	}

	/**
	 * Check 
	 * @param str The string to check
	 * @returns Whether true of the string has a valid length
	 */
	Boolean IsValidLineLength(String str)
	{
		return str.length() == LineLength() ? true : false;
	}

	/**
	 * Convert a string containing an integer
	 * @param[in] str The string to convert
	 * @param[out] val The result
	 * @exception TleException on conversion error
	 */
	int ExtractInteger(String str)
	{
		return Integer.parseInt(str.trim());
	}

	/**
	 * Convert a string containing an double
	 * @param[in] str The string to convert
	 * @param[in] point_pos The position of the decimal point. (-1 if none)
	 * @param[out] val The result
	 * @exception TleException on conversion error
	 */
	double ExtractDouble(String str, int point_pos)
	{
		String tmp = new String(str);


		if (point_pos==-1 && (tmp.charAt(0)=='-'|| tmp.charAt(0)=='+'))
		{
			tmp = tmp.charAt(0)+"."+tmp.substring(1);
		}
		else if (point_pos==-1)
		{
			tmp= '.'+tmp;
		}
		return  Double.parseDouble(tmp);
	}

	/**
	 * Convert a string containing an exponential
	 * @param[in] str The string to convert
	 * @param[out] val The result
	 * @exception TleException on conversion error
	 */
	double  ExtractExponential(String str) throws TleException
	{
		String tmp = new String();
		for (int i =0;i<str.length();i++)
		{
			if (i==0)
			{
				if ((str.charAt(i) == '-' || str.charAt(i) == '+'||str.charAt(i) == ' '))
				{
					if (str.charAt(i)=='-')
						tmp+='-';
					tmp+='0';
					tmp+='.';
				}
				else
				{
					throw new TleException("Invalid sign");
				}
			}
			else if (i==str.length()-2)
			{
				if (str.charAt(i) == '-' || str.charAt(i)=='+')
				{
					tmp+='e';
					tmp+=str.charAt(i);
				}
			}
			else
			{
	            if (Character.isDigit(str.charAt(i)))
	            {
	                tmp += str.charAt(i);
	            }
	            else
	            {
	                throw new TleException("Invalid digit");
	            }
			}
		}


		return  Double.parseDouble(tmp);
	}
}