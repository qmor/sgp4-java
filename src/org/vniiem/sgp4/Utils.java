package org.vniiem.sgp4;

import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.temporal.ChronoUnit;
import java.util.Calendar;

public class Utils {
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



	/*
	 * always positive result
	 * Mod(-3,4)= 1   fmod(-3,4)= -3
	 */
	public static double Mod( double x,  double y)
	{
		if (y == 0.0)
		{
			return x;
		}

		return x - y * Math.floor(x / y);
	}

	public static double WrapNegPosPI(double a)
	{
		return Mod(a + Globals.kPI, Globals.kTWOPI) - Globals.kPI;
	}

	public static double WrapTwoPI( double a)
	{
		return Mod(a, Globals.kTWOPI);
	}

	double WrapNegPos180( double a)
	{
		return Mod(a + 180.0, 360.0) - 180.0;
	}

	public static double TotalMinutes(LocalDateTime dt)
	{
		return Utils.getDoubleUnixTime(dt)/60.0;
	}

	public static double Wrap360( double a)
	{
		return Mod(a, 360.0);
	}

	static double DegreesToRadians( double degrees)
	{
		return degrees * Globals.kPI / 180.0;
	}

	double RadiansToDegrees( double radians)
	{
		return radians * 180.0 / Globals.kPI;
	}

	public static double AcTan( double sinx,  double cosx)
	{
		if (cosx == 0.0)
		{
			if (sinx > 0.0)
			{
				return Globals.kPI / 2.0;
			}
			else
			{
				return 3.0 * Globals.kPI / 2.0;
			}
		}
		else
		{
			if (cosx > 0.0)
			{
				return Math.atan(sinx / cosx);
			}
			else
			{
				return Globals.kPI + Math.atan(sinx / cosx);
			}
		}
	}

	public static double ToLocalMeanSiderealTime(LocalDateTime dt,final double lon)
	{
		return Utils.WrapTwoPI(ToGreenwichSiderealTime(dt) + lon);
	}

	public static double ToGreenwichSiderealTime(LocalDateTime dt) 
	{

		// t = Julian centuries from 2000 Jan. 1 12h UT1
		//todo: checkit
		final double t = (Utils.toJulianDay(Utils.getDoubleUnixTime(dt)) - 2451545.0) / 36525.0;

		// Rotation angle in arcseconds
		double theta = 67310.54841
				+ (876600.0 * 3600.0 + 8640184.812866) * t
				+ 0.093104 * t * t
				- 0.0000062 * t * t * t;

		// 360.0 / 86400.0 = 1.0 / 240.0
		return Utils.WrapTwoPI(Utils.DegreesToRadians(theta / 240.0));
	}

	public static LocalDateTime epochDateTime(int year, double day)
	{
		LocalDateTime dt =  LocalDateTime.of(year,1,1,0,0);
		dt = dt.plus((long)((day-1)*86400.0*1e9), ChronoUnit.NANOS);
		return dt;
	}
	public static String SubString(String input, int startindex, int length)
	{
		return new String(input.getBytes(),startindex,length);
	}

	public static double getDoubleUnixTime(LocalDateTime dt)
	{
		Instant helper =dt.toInstant(ZoneOffset.ofTotalSeconds(0)); 
		return helper.getEpochSecond()+helper.getNano()/1e9;
	}

	public static double toJulianDay(double unixtime)
	{
		return ( unixtime / 86400.0 ) + 2440587.5;
	}
	
	public static LocalDateTime minusdt(LocalDateTime d1, LocalDateTime d2)
	{
		double dt1 = Utils.getDoubleUnixTime(d1);
		double dt2 = Utils.getDoubleUnixTime(d2);
		
		double res = dt1-dt2;
		LocalDateTime result = LocalDateTime.of(1970, 1, 1, 0, 0 );
		result = result.plus((long)(res*1e9), ChronoUnit.NANOS);
		return result;
		
	}

	public static LocalDateTime plusdt(LocalDateTime d1, LocalDateTime d2)
	{
		double dt1 = Utils.getDoubleUnixTime(d1);
		double dt2 = Utils.getDoubleUnixTime(d2);
		
		double res = dt1+dt2;
		LocalDateTime result = LocalDateTime.of(1970, 1, 1, 0, 0 );
		result = result.plus((long)(res*1e9), ChronoUnit.NANOS);
		return result;
		
	}

}
