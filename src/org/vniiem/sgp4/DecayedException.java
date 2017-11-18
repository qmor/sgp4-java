package org.vniiem.sgp4;

import java.time.LocalDateTime;

import javax.vecmath.Vector4d;

public class DecayedException extends Exception{

	
	  public DecayedException( LocalDateTime dt,  Vector4d pos,  Vector4d vel)
      
  {
		  super("Satellite decayed");
		   _dt =(dt);
	       _pos=(pos);
	       _vel=(vel);
  }

  /**
   * @returns the date
   */
  LocalDateTime Decayed() 
  {
      return _dt;
  }

  /**
   * @returns the position
   */
  Vector4d Position() 
  {
      return _pos;
  }

  /**
   * @returns the velocity
   */
  Vector4d Velocity() 
  {
      return _vel;
  }


  LocalDateTime _dt;
  Vector4d _pos;
  Vector4d _vel;
	
	
}
