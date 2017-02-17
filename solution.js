//
// Working solution to http://hellospace.reaktor.com/
//
// Antti Kangasrääsiö
// antti.kangasraasio@iki.fi
// 17.2.2017
//

function getPlanetState(state, planetName) {
  var planet;
  for (planet in state.planetStates) {
    var ps = state.planetStates[planet];
    if (ps.name == planetName) {
      return ps
    }
  }
  return null;
}

function pointingToDir(target, object) {
  // x-axis is forward
  var dirTarget = object.rotation.inverse().vmult(target).unit();
  return dirTarget.dot(new Vec3(1,0,0)) > 0.95;
}

function changeApoapsis(targetRadius, rocket, planet) {
  var orbit = getOrbit(rocket, planet);
  var precision = 2.0;
  var thrust = 0.0;
  var targetDir = orbit.relVelocity;
  var msg = "";
  if (orbit.R_a < targetRadius - precision) {
    msg = "Attempting to raise apoapsis to " + targetRadius;
    if (Math.abs(orbit.r - orbit.R_p) < precision) {
      msg = msg + ", at periapsis";
      if (pointingToDir(targetDir, rocket)) {
        msg = msg + ", thrusting!";
        thrust = 1.0;
      }
    }
  }
  else if (orbit.R_a > targetRadius + precision) {
    msg = "Attempting to lower apoapsis to " + targetRadius;
    targetDir = targetDir.clone().negate();
    if (Math.abs(orbit.r - orbit.R_p) < precision) {
      msg = msg + ", at periapsis";
      if (pointingToDir(targetDir, rocket)) {
        msg = msg + ", thrusting!";
        thrust = 1.0;
      }
    }
  }
  else {
    msg = "Apoapsis OK";
  }
  console.log(msg);
  var rcs = controlDir(targetDir, rocket);
  return {
      thrust: thrust,
      rcs: rcs}
}

function goToward(dir, rocket) {
  var rcs = controlDir(dir, rocket);
  return {
      thrust: 0.0,
      rcs: rcs}
}

function changePeriapsis(targetRadius, rocket, planet) {
  var orbit = getOrbit(rocket, planet);
  var precision = 2.0;
  var thrust = 0.0;
  var targetDir = orbit.relVelocity;
  var msg = "";
  if (orbit.R_p < targetRadius - precision) {
    msg = "Attempting to raise periapsis to " + targetRadius;
    if (Math.abs(orbit.r - orbit.R_a) < precision) {
      msg = msg + ", at apoapsis";
      if (pointingToDir(targetDir, rocket)) {
        console.log("Pointing forward, thrusting");
        thrust = 1.0;
      }
    }
  }
  else if (orbit.R_p > targetRadius + precision) {
    msg = "Attempting to lower periapsis to " + targetRadius;
    targetDir = targetDir.clone().negate();
    if (Math.abs(orbit.r - orbit.R_a) < precision) {
      msg = msg + ", at apoapsis";
      if (pointingToDir(targetDir, rocket)) {
        msg = msg + ", thrusting!";
        thrust = 1.0;
      }
    }
  }
  else {
    msg = "Apoapsis OK";
  }
  console.log(msg);
  var rcs = controlDir(targetDir, rocket);
  return {
      thrust: thrust,
      rcs: rcs}
}

function roundOrbit(orbit, targetOrbitRadius, tolerance, rocket, planet) {
  if (Math.abs(orbit.R_a - targetOrbitRadius) > tolerance) {
    return changeApoapsis(targetOrbitRadius, rocket, planet);
  }
  else if (Math.abs(orbit.R_p - targetOrbitRadius) > tolerance) {
    return changePeriapsis(targetOrbitRadius, rocket, planet);
  }
  return goToward(orbit.relVelocity, rocket);
}

function launch(targetDir, rocket, planet) {
  var rcs = controlDir(targetDir, rocket);
  return {
      thrust: 1.0,
      rcs: rcs}
}

function land(targetDir, orbit, maxVelocity, rocket, planet) {
  var rcs = controlDir(targetDir, rocket);
  var thrust = 0.0;
  if (orbit.v > maxVelocity) {
    thrust = 1.0;
  }
  return {
      thrust: thrust,
      rcs: rcs}
}

function controlDir(target, rocket) {
  // target in absolute coords
  // orientation control
  var dirTarget = rocket.rotation.inverse().vmult(target).unit();
  var p_o = 2.0;
  var pitch = p_o*dirTarget.z;
  var yaw = p_o*dirTarget.y;
  var roll = 0.0;
  // angular velocity stabilizer
  var p_v = 1.0;
  pitch = pitch - p_v*rocket.angularVelocity.y;
  yaw = yaw - p_v*rocket.angularVelocity.z;
  roll = roll - p_v*rocket.angularVelocity.x;
  // rate limiter
  return {pitch: Math.max(-0.5, Math.min(0.5, pitch)),
           yaw: Math.max(-0.5, Math.min(0.5, yaw)),
           roll: Math.max(-0.5, Math.min(0.5, roll))}
}

function getG() {
  return 3.2173623151255695e-9;
  // they did not document the world params..
  var moonApoapsisRadius = 298.908;
  var moonApoapsisVelocity = 10.053;
  var moonPeriapsisRadius = 280.155;
  var moonPeriapsisVelocity = 10.763;
  var moonMass = 730000000000;
  var earthMass = 5900000000000;
  var a = moonApoapsisRadius + moonPeriapsisRadius;
  var e = 1 - moonPeriapsisRadius / a
  var h1 = moonApoapsisRadius*moonApoapsisVelocity;
  var h2 = moonPeriapsisRadius*moonPeriapsisVelocity;
  var h = (h1+h2)/2.0;
  var mu = Math.pow(h, 2) / (a * (1 - Math.pow(e, 2)));
  var G = mu / (moonMass + earthMass);
  return G
}

function getOrbit(orbiter, planet) {
  var relPosition = orbiter.position.vsub(planet.position);
  var posDir = relPosition.unit();
  var relVelocity = orbiter.velocity.vsub(planet.velocity);
  var velDir = relVelocity.unit();
  var gamma = Math.acos(posDir.dot(velDir));
  var G = getG();
  var GM = G * (planet.mass + orbiter.mass);
  var r = orbiter.position.distanceTo(planet.position);
  var v = relVelocity.length();
  var v2 = relVelocity.lengthSquared()
  var e = Math.sqrt(
    Math.pow(r * v2 / GM - 1, 2)
      * Math.pow(Math.sin(gamma), 2)
      + Math.pow(Math.cos(gamma), 2));
  var a = 1.0 / ( 2.0 / r - v2 / GM );
  var R_p = a * (1 - e);
  var R_a = a * (1 + e);
  return {
          r: r,
          v: v,
          a: a,
          e: e,
          R_p: R_p,
          R_a: R_a,
          relPosition: relPosition,
          relVelocity: relVelocity,
         }
}

return function GoToMoon(state) {
 var moon = getPlanetState(state, "Moon");
 var earth = getPlanetState(state, "Earth2");
 var rocket = state.rocket;
   
 var orbitE = getOrbit(rocket, earth);
 var orbitM = getOrbit(rocket, moon);

 if (orbitM.r > 50 || orbitM.v > 10) {
   console.log("E: r",orbitE.r,"v",orbitE.v,"Rp",orbitE.R_p,"Ra",orbitE.R_a)
   var moonOrbitRadius = moon.position.distanceTo(earth.position);
   var rocketToMoon = moon.position.vsub(rocket.position).unit();
     
   if (orbitE.r < 30.0) {
     console.log("Launch")
     return new Controls(launch(rocketToMoon, rocket, earth));
   }
   else if (orbitE.R_p < 180) {
     console.log("Low earth orbit")
     return new Controls(changePeriapsis(200, rocket, earth));
   }
   else {
     console.log("To Moon orbit")
     return new Controls(roundOrbit(orbitE, moonOrbitRadius, 40, rocket, earth));
   }
 } else {
   console.log("M: r",orbitM.r,"v",orbitM.v,"Rp",orbitM.R_p,"Ra",orbitM.R_a);
   var moonToRocket = rocket.position.vsub(moon.position).unit();
    
   if (orbitM.v > 5) {
     console.log("High Moon orbit")
     var targetOrbitRadius = moon.radius + 10;
     return new Controls(roundOrbit(orbitM, targetOrbitRadius, 3, rocket, moon));
   } else {
     console.log("Land to Moon")
     return new Controls(land(moonToRocket, orbitM, 2, rocket, moon));
   }
 }
}
