package org.vniiem.sgp4;



public class CoordGeodetic {


	public CoordGeodetic()
	{
		
	}

	    /**
	     * Constructor
	     * @param[in] lat the latitude (degrees by default)
	     * @param[in] lon the longitude (degrees by default)
	     * @param[in] alt the altitude in kilometers
	     * @param[in] is_radians whether the latitude/longitude is in radians
	     */
	    public CoordGeodetic(
	            double lat,
	            double lon,
	            double alt,
	            Boolean is_radians)
	    {
	        if (is_radians)
	        {
	            latitude = lat;
	            longitude = lon;
	        }
	        else
	        {
	            latitude = Utils.DegreesToRadians(lat);
	            longitude = Utils.DegreesToRadians(lon);
	        }
	        altitude = alt;
	    }

	    /**
	     * Copy constructor
	     * @param[in] geo object to copy from
	     */
	   public CoordGeodetic(CoordGeodetic geo)
	    {
	        latitude = geo.latitude;
	        longitude = geo.longitude;
	        altitude = geo.altitude;
	    }


	    /**
	     * Dump this object to a string
	     * @returns string
	     */
	    @Override
		public   String toString() 
	    {
	    	/*
	        std::stringstream ss;
	        ss << std::right << std::fixed << std::setprecision(3);
	        ss << "Lat: " << std::setw(7) << Util::RadiansToDegrees(latitude);
	        ss << ", Lon: " << std::setw(7) << Util::RadiansToDegrees(longitude);
	        ss << ", Alt: " << std::setw(9) << altitude;
	        return ss.str();
	        */
	        return super.toString();
	    }

	    /** latitude in radians (-PI >= latitude < PI) */
	   public double latitude;
	    /** latitude in radians (-PI/2 >= latitude <= PI/2) */
	   public double longitude;
	    /** altitude in kilometers */
	   public double altitude;
	}



