package org.vniiem.sgp4;

import java.time.LocalDateTime;

import javax.vecmath.Vector3d;
import javax.vecmath.Vector4d;

public class Eci {

	/**
	 * @param[in] dt the date to be used for this position
	 * @param[in] latitude the latitude in degrees
	 * @param[in] longitude the longitude in degrees
	 * @param[in] altitude the altitude in kilometers
	 */
	public Eci( LocalDateTime dt,
			final double latitude,
			final double longitude,
			final double altitude)
	{
		ToEci(dt, new CoordGeodetic(latitude, longitude, altitude,false));
	}

	/**
	 * @param[in] dt the date to be used for this position
	 * @param[in] geo the position
	 */
	Eci(final LocalDateTime dt, final CoordGeodetic geo)
	{
		ToEci(dt, geo);
	}

	/**
	 * @param[in] dt the date to be used for this position
	 * @param[in] position the position
	 */
	Eci(final LocalDateTime dt, final Vector4d position)


	{
		m_dt = (dt);
		m_position = (position);
	}


	/**
	 * @param[in] dt the date to be used for this position
	 * @param[in] position the position
	 * @param[in] velocity the velocity
	 */
	public Eci(final LocalDateTime dt, final Vector4d position, final Vector4d velocity)

	{
		m_dt =(dt);
		m_position = (position);
		m_velocity = (velocity);

	}



	/**
	 * Update this object with a new date and geodetic position
	 * @param dt new date
	 * @param geo new geodetic position
	 */
	void Update(final LocalDateTime dt, final CoordGeodetic geo)
	{
		ToEci(dt, geo);
	}

	/**
	 * @returns the position
	 */
	public Vector4d Position()
	{
		return m_position;
	}

	/**
	 * @returns the velocity
	 */
	public Vector4d Velocity()
	{
		return m_velocity;
	}
	
	/**
	 * @returns the position in array of doubles
	 */
	public double[] PositionAsDoubleArray()
	{
		return new double[] {m_position.x, m_position.y, m_position.z};
	}

	/**
	 * @returns the velocity in array of doubles
	 */
	public double[] VelocityAsDoubleArray()
	{
		return new double[] {m_velocity.x, m_velocity.y, m_velocity.z};
	}

	/**
	 * @returns the date
	 */
	LocalDateTime GetDateTime()
	{
		return m_dt;
	}

	/**
	 * @returns the position in geodetic form
	 */



	LocalDateTime m_dt;

	Vector4d m_position;
	Vector4d m_velocity;
	static final double mfactor = Globals.kTWOPI * (Globals.kOMEGA_E / Globals.kSECONDS_PER_DAY);
	void ToEci(final LocalDateTime dt, final CoordGeodetic geo)
	{
		/*
		 * set date
		 */
		m_dt = dt;
		/*
		 * Calculate Local Mean Sidereal Time for observers longitude
		 */

		final double theta = Utils.ToLocalMeanSiderealTime(m_dt,geo.longitude);

		/*
		 * take into account earth flattening
		 */
		final double c = 1.0
				/ Math.sqrt(1.0 + Globals.kF * (Globals.kF - 2.0) * Math.pow(Math.sin(geo.latitude), 2.0));
		final double s = Math.pow(1.0 - Globals.kF, 2.0) * c;
		final double achcp = (Globals.kXKMPER * c + geo.altitude) * Math.cos(geo.latitude);

		/*
		 * X position in km
		 * Y position in km
		 * Z position in km
		 * W magnitude in km
		 */
		m_position.x = achcp * Math.cos(theta);
		m_position.y = achcp * Math.sin(theta);
		m_position.z = (Globals.kXKMPER * s + geo.altitude) * Math.sin(geo.latitude);
		m_position.w = m_position.length();

		/*
		 * X velocity in km/s
		 * Y velocity in km/s
		 * Z velocity in km/s
		 * W magnitude in km/s
		 */
		m_velocity.x = -mfactor * m_position.y;
		m_velocity.y = mfactor * m_position.x;
		m_velocity.z = 0.0;
		m_velocity.w = m_velocity.length();
	}

	/**
	 * @returns the position in geodetic form
	 */  
	static final double e2 = Globals.kF * (2.0 - Globals.kF);

	CoordGeodetic ToGeodetic()
	{
		final double theta = Utils.AcTan(m_position.y, m_position.x);

		final double lon = Utils.WrapNegPosPI(theta - Utils.ToGreenwichSiderealTime(m_dt));

		final double r = Math.sqrt((m_position.x * m_position.x) + (m_position.y * m_position.y));


		double lat = Utils.AcTan(m_position.z, r);
		double phi = 0.0;
		double c = 0.0;
		int cnt = 0;

		do
		{
			phi = lat;
			final double sinphi = Math.sin(phi);
			c = 1.0 / Math.sqrt(1.0 - e2 * sinphi * sinphi);
			lat = Utils.AcTan(m_position.z + Globals.kXKMPER * c * e2 * sinphi, r);
			cnt++;
		}
		while (Math.abs(lat - phi) >= 1e-10 && cnt < 10);

		final double alt = r / Math.cos(lat) - Globals.kXKMPER * c;

		return new CoordGeodetic(lat, lon, alt, true);
	}
}