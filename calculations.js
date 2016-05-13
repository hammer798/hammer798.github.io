// JavaScript Document
function getData(id){
	id = id.toString();
	id = document.getElementById(id).value;
	id = parseFloat(id);
	return id;
}

function rotKinematics1(){
	var deltatheta = getData('deltatheta');
	var initOmega = getData('initomega');
	var alpha = getData('alpha');
	var time = getData('rottime');
	if(isNaN(deltatheta)){
		deltatheta = initOmega*time + .5*alpha*Math.pow(time, 2);
		deltatheta = "\u0394\u03B8 = " + deltatheta.toPrecision(3) + " rads";
		return deltatheta;
	}
	else if(isNaN(initOmega)){
		initOmega = (deltatheta - .5*alpha*Math.pow(time, 2))/time;
		initOmega = "\u03C90 = " + initOmega.toPrecision(3) + " rad/s";
		return initOmega;
	}
	else if(isNaN(alpha)){
		alpha = 2*(deltatheta - initOmega*time)/Math.pow(time, 2);
		alpha = "\u03B1 = " + alpha.toPrecision(3) + " rad/s^2";
		return alpha;
	}
	else{
		if(initOmega === 0){
			time = Math.sqrt(2*deltatheta/alpha);
		}
		else{
			time = (-1*initOmega + Math.sqrt(Math.pow(initOmega, 2) + 2*alpha*deltatheta))/alpha;
		}
		time = "t = " + time.toPrecision(3) + " s";
		return time;
	}

}
function rotKinematics2(){
	var omega = getData('omega');
	var initomega = getData('initomega');
	var time = getData('rottime');
	var alpha = getData('alpha');
	if(isNaN(omega)){
		omega = initomega + alpha*time;
		omega = "\u03C9 = " + omega.toPrecision(3) + " rad/s";
		return omega;
	}
	else if(isNaN(initomega)){
		initomega = omega - alpha*time;
		initomega = "\u03C90 = " + initomega.toPrecision(3) + " rad/s";
		return initomega;
	}
	else if(isNaN(alpha)){
		alpha = (omega - initomega)/time;
		alpha = "\u03B1 = " + alpha.toPrecision(3) + " rad/s^2";
		return alpha;
	}
	else{
		time = (omega - initomega)/alpha;
		time = "t = " + time.toPrecision(3) + " s";
		return time;
	}
}
function rotKinematics3(){
		var centAccel = getData('centaccel');
		var rotvel = getData('rotvel');
		var radius = getData('radius');
		if(isNaN(centAccel)){
			centAccel = Math.pow(rotvel, 2)/radius;
			centAccel = "ac = " + centAccel.toPrecision(3) + " m/s^2";
			return centAccel;
		}
		else if(isNaN(rotvel)){
			rotvel = Math.sqrt(centAccel*radius);
			rotvel = "v = " + rotvel.toPrecision(3) + " m/s";
			return rotvel;
		}
		else{
			radius = Math.pow(rotvel, 2)/centAccel;
			radius = "r = " + radius.toPrecision(3) + " m";
			return radius;
		}
}

function chooseRotFormula(){
	var which = document.getElementById("select2");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "rotone"){
		answer = rotKinematics1();
	}
	else if(formula === "rottwo"){
		answer = rotKinematics2();	
	}
	else{
		answer = rotKinematics3();
	}
	document.getElementById("rotanswer").innerHTML = answer;
}
function changeRotVisible(){
	var which = document.getElementById("select2");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["omegapara", "thetapara", "alphapara", "rottimepara", "initomegapara", "rotvelpara", "centaccelpara", "radpara"];
	var ids2 = ["deltatheta", "initomega", "omega", "alpha", "rottime", "rotvel", "centaccel", "radius"];
	var settings = ["none", "inline", "inline", "inline", "inline", "none", "none", "none"];
	var settings2 = ["inline", "none", "inline", "inline", "inline", "none", "none", "none"];
	var settings3 = ["none", "none", "none", "none", "none", "inline", "inline" ,"inline"];
	if(formula === "rotone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];
		}
	}
	else if(formula === "rottwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];	
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";	
	}
	document.getElementById("rotanswer").innerHTML = "";
}
function kinematics1(){
	var deltax = getData('deltax');
	var v0 = getData('initvelocity');
	var a = getData('acceleration');
	var t = getData('time');
	if(isNaN(deltax)){
		deltax = v0*t + (.5*a*Math.pow(t, 2));
		deltax = "\u0394x = " + deltax.toPrecision(3) + " m";
		return deltax;			
	}
	else if(isNaN(v0)){
		v0 = (deltax - (.5*a*Math.pow(t, 2)))/t;
		v0 = "v0 = " + v0.toPrecision(3) + " m/s";	
		return v0;
	}
	else if(isNaN(a)){
		a = (2*(deltax - (v0*t)))/Math.pow(t, 2);
		a = "a = " + a.toPrecision(3) + " m/s^2";
		return a;	
	}
	else{
		if(v0 === 0){
			t = Math.sqrt(((2*deltax)/a));	
		}
		else{
			t = (v0*-1 + Math.sqrt(Math.pow(v0, 2) - 2*a*deltax*-1))/a;	
			
		}
		t = "t = " + t.toPrecision(3) + " s";
		return t;
	}
	
}
function kinematics2(){
	var velocity = getData('velocity');
	var initvelocity = getData('initvelocity');
	var acceleration = getData('acceleration');
	var time = getData('time');
	if(isNaN(velocity)){
		velocity = initvelocity + (acceleration*time);
		velocity = "v = " + velocity.toPrecision(3) + " m/s";
		return velocity;
	}
	else if(isNaN(initvelocity)){
		initvelocity = velocity - acceleration*time;
		initvelocity = "v0 = " + initvelocity.toPrecision(3) + " m/s";
		return initvelocity;	
	}
	else if(isNaN(acceleration)){
		acceleration = (velocity - initvelocity)/time;
		acceleration = "a = " + acceleration.toPrecision(3) + " m/s^2";
		return acceleration;
	}
	else{
		time = (velocity - initvelocity)/acceleration;
		time = "t = " + time.toPrecision(3) + " s";
		return time;	
	}
}
function kinematics3(){
	var velocity = getData('velocity');
	var initvelocity = getData('initvelocity');
	var deltax = getData('deltax');
	var acceleration = getData('acceleration');
	if(isNaN(velocity)){
		velocity = Math.sqrt(Math.pow(initvelocity, 2) + 2*acceleration*deltax);
		velocity = "v = " + velocity.toPrecision(3) + " m/s";
		return velocity;	
	}
	else if(isNaN(initvelocity)){
		initvelocity = Math.sqrt(Math.pow(velocity, 2) - 2*acceleration*deltax);
		initvelocity = "v0 = " + initvelocity.toPrecision(3) + " m/s";
		return initvelocity;	
	}
	else if(isNaN(deltax)){
		deltax = (Math.pow(velocity, 2) - Math.pow(initvelocity, 2))/(2*acceleration);
		deltax = "\u0394x = " + deltax.toPrecision(3) + " m";
		return deltax;
	}
	else{
		acceleration = (Math.pow(velocity, 2) - Math.pow(initvelocity, 2))/(2*deltax);
		acceleration = "a = " + acceleration.toPrecision(3) + " m/s^2";
		return acceleration;	
	}
}
function chooseFormula(){
	var which = document.getElementById("select");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "kinone"){
		answer = kinematics1();
	}
	else if(formula === "kintwo"){
		answer = kinematics2();	
	}
	else{
		answer = kinematics3();	
	}
	document.getElementById("answer").innerHTML = answer;
}
function changeVisible(){
	var which = document.getElementById("select");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["velocity", "initvelocity", "acceleration", "deltax", "time"];
	if(formula === "kinone"){
		document.getElementById("velpara").style.display = "none";
		document.getElementById("delpara").style.display = "inline";
		document.getElementById("timepara").style.display = "inline";
	}	
	else if(formula === "kintwo"){
		document.getElementById("delpara").style.display = "none";
		document.getElementById("timepara").style.display = "inline";
		document.getElementById("velpara").style.display = "inline";
			
	}
	else{
		document.getElementById("timepara").style.display = "none";	
		document.getElementById("velpara").style.display = "inline";
		document.getElementById("delpara").style.display = "inline";
	}
	for(var i in ids){
		document.getElementById(ids[i]).value = "";	
	}
	document.getElementById("answer").innerHTML = "";
}
function netForce(){
	var fnet = getData('netforce');
	var mass = getData('mass');
	var acceleration = getData('foraccel');
	if(isNaN(fnet)){
		fnet = mass*acceleration;
		fnet = "Fnet = " + fnet.toPrecision(3) + " N";
		return fnet;
	}
	else if(isNaN(mass)){
		mass = fnet/acceleration;
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else{
		acceleration = fnet/mass;
		acceleration = "a = " + acceleration.toPrecision(3) + " m/s^2";
		return acceleration;
	}
}
function forceFriction(){
	var forceFric = getData('forcefric');
	var mu = getData('mu');
	var forceNorm = getData('normforce');
	if(isNaN(forceFric)){
		forceFric = mu*forceNorm;
		forceFric = "Ff = " + forceFric.toPrecision(3) + " N";
		return forceFric;
	}
	else if(isNaN(mu)){
		mu = forceFric/forceNorm;
		mu = "\u03BC = " + mu.toPrecision(3);
		return mu;
	}
	else{
		forceNorm = forceFric/mu;
		forceNorm = "Fn = " + forceNorm.toPrecision(3) + " N";
		return forceNorm;
	}
}
function forceGrav(){
	var gravForce = getData('gravforce');
	var mass = getData('mass');
	if(isNaN(gravForce)){
		gravForce = mass*9.8;
		gravForce = "Fg = " + gravForce.toPrecision(3) + " N";
		return gravForce;
	}
	else{
		mass = gravForce/9.8;
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
}
function univGrav(){
	var gravForce = getData('gravforce');
	var mass = getData('mass');
	var masstwo = getData('masstwo');
	var radius = getData('forcerad');
	if(isNaN(gravForce)){
		gravForce = (6.67e-11*mass*masstwo)/Math.pow(radius, 2);
		gravForce = "Fg = " + gravForce.toPrecision(3) + " N";
		return gravForce;
	}
	else if(isNaN(mass)){
		mass = (gravForce*Math.pow(radius, 2))/(6.67e-11*masstwo);
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else if(isNaN(masstwo)){
		masstwo = (gravForce*Math.pow(radius, 2))/(6.67e-11*mass);
		masstwo = "m = " + masstwo.toPrecision(3) + " kg";
		return masstwo;
	}
	else{
		radius = Math.sqrt((6.67e-11*mass*masstwo)/gravForce);
		radius = "r = " + radius.toPrecision(3) + " m";
		return radius;
	}
}
function springForce(){
	var forceSpring = getData('forcespring');
	var k = getData('k');
	var stretch = getData('xstretch');
	if(isNaN(forceSpring)){
		forceSpring = k*stretch;
		forceSpring = "Fs = " + forceSpring.toPrecision(3) + " N";
		return forceSpring;
	}
	else if(isNaN(k)){
		k = forceSpring/stretch;
		k = "k = " + k.toPrecision(3) + " N/m";
		return k;
	}
	else{
		stretch = forceSpring/k;
		stretch = "\u0394x = " + stretch.toPrecision(3) + " m";
		return stretch;
	}
}
function chooseForceFormula(){
	var which = document.getElementById("forceselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "forceone"){
		answer = netForce();
	}
	else if(formula === "forcetwo"){
		answer = forceFriction();	
	}
	else if(formula === "forcethree"){
		answer = forceGrav();	
	}
	else if(formula === "forcefour"){
		answer = univGrav();
	}
	else{
		answer = springForce();
	}
	document.getElementById("foranswer").innerHTML = answer;
}
function changeForceVisible(){
	var which = document.getElementById("forceselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["netforpara", "masspara", "foraccelpara", "forcefricpara", "mupara", "normforcepara", "gravforcepara", "masstwopara", "forceradpara", "forcesprpara", "kpara", "xstretchpara"];
	var ids2 = ["netforce", "mass", "foraccel", "forcefric", "mu", "normforce", "gravforce", "masstwo", "forcerad", "forcespring", "k", "xstretch"];
	var settings = ["inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none"];
	var settings3 = ["none", "inline", "none", "none", "none", "none", "inline", "none", "none", "none", "none", "none"];
	var settings4 = ["none", "inline", "none", "none", "none", "none", "inline", "inline", "inline", "none", "none", "none"];
	var settings5 = ["none", "none", "none", "none", "none", "none", "none", "none", "none", "inline", "inline", "inline"];
	if(formula === "forceone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];	
		}
	}	
	else if(formula === "forcetwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];	
		}
	}
	else if(formula === "forcethree"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];	
		}
	}
	else if(formula === "forcefour"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings4[i];	
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings5[i];	
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";
	}
	document.getElementById("answer").innerHTML = "";
}
function momentum1(){
	var momentum = getData('moment');
	var mass = getData('mommass');
	var velocity = getData('momvel');
	if(isNaN(momentum)){
		momentum = mass*velocity;
		momentum = "p = " + momentum.toPrecision(3) + " kgm/s";
		return momentum;
	}
	else if(isNaN(mass)){
		mass = momentum/velocity;
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else{
		velocity = momentum/mass;
		velocity = "v = " + velocity.toPrecision(3) + " m/s";
		return velocity;
	}
}
function momentum2(){
	var impulse = getData('deltamoment');
	var forceNet = getData('momforce');
	var time = getData('momtime');
	if(isNaN(impulse)){
		impulse = forceNet*time;
		impulse = "\u0394p = " + impulse.toPrecision(3) + " kgm/s";
		return impulse;
	}
	else if(isNaN(forceNet)){
		forceNet = impulse/time;
		forceNet = "Fnet = " + forceNet.toPrecision(3) + " N";
		return forceNet;
	}
	else{
		time = impulse/forceNet;
		time = "t = " + time.toPrecision(3) + " s";
		return time;
	}
}
function chooseMomFormula(){
	var which = document.getElementById("momentselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "momentumone"){
		answer = momentum1();
	}
	else{
		answer = momentum2();	
	}
	document.getElementById("momanswer").innerHTML = answer;
}
function changeMomVisible(){
	var which = document.getElementById("momentselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["momentpara", "mommasspara", "momvelpara", "deltamomentpara", "momforcepara", "momtimepara"];
	var ids2 = ["moment", "mommass", "momvel", "deltamoment", "momforce", "momtime"];
	var settings = ["inline", "inline", "inline", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "inline", "inline", "inline"];
	if(formula === "momentumone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];	
		}
	}	
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];	
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";	
	}
	document.getElementById("momanswer").innerHTML = "";
}
function kineticEnergy(){
	var kineticEnergy = getData('kinener');
	var mass = getData('enermass');
	var velocity = getData('enervel');
	if(isNaN(kineticEnergy)){
		kineticEnergy = .5*mass*Math.pow(velocity, 2);
		kineticEnergy = "KE = " + kineticEnergy.toPrecision(3) + " J";
		return kineticEnergy;
	}
	else if(isNaN(mass)){
		mass = kineticEnergy*2/Math.pow(velocity, 2);
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else{
		velocity = Math.sqrt(kineticEnergy*2/mass);
		velocity = "v = " + velocity.toPrecision(3) + " m/s";
		return velocity;
	}
}

function calcWork(){
	var work = getData('work');
	var force = getData('workforce');
	var distance = getData('distance');
	if(isNaN(work)){
		work = force*distance;
		work = "W = " + work.toPrecision(3) + " J";
		return work;
	}
	else if(isNaN(force)){
		force = work/distance;
		force = "F = " + force.toPrecision(3) + " N";
		return force;
	}
	else{
		distance = work/force;
		distance = "d = " + distance.toPrecision(3) + " m";
		return distance;
	}
}

function calcPower(){
	var power = getData('power');
	var deltaEnergy = getData('enerchange');
	var time = getData('enertime');
	if(isNaN(power)){
		power = deltaEnergy/time;
		power = "P = " + power.toPrecision(3) + " W";
		return power;
	}
	else if(isNaN(deltaEnergy)){
		deltaEnergy = power*time;
		deltaEnergy = "\u0394E = " + deltaEnergy.toPrecision(3) + " J";
		return deltaEnergy;
	}
	else{
		time = deltaEnergy/power;
		time = "t = " + time.toPrecision(3) + " s";
		return time;
	}
}
function rotKinetic(){
	var rotKineticEn = getData('rotkin');
	var rotInert = getData('rotinert');
	var omega = getData('eneromega');
	if(isNaN(rotKineticEn)){
		rotKineticEn = .5*rotInert*Math.pow(omega, 2);
		rotKineticEn = "KE = " + rotKineticEn.toPrecision(3) + " J";
		return rotKineticEn;
	}
	else if(isNaN(rotInert)){
		rotInert = rotKineticEn*2/Math.pow(omega, 2);
		rotInert = "I = " + rotInert.toPrecision(3) + " kgm^2";
		return rotInert;
	}
	else{
		omega = Math.sqrt(rotKineticEn*2/rotInert);
		omega = "\u03C9	= " + omega.toPrecision(3) + " rad/s";
		return omega;
	}
}
function energySpring(){
	var springEnergy = getData('enerspr');
	var k = getData('enerk');
	var x = getData('enerstretch');
	if(isNaN(springEnergy)){
		springEnergy = .5*k*Math.pow(x, 2);
		springEnergy = "Us = " + springEnergy.toPrecision(3) + " J";
		return springEnergy;
	}
	else if(isNaN(k)){
		k = springEnergy*2/Math.pow(x, 2);
		k = "k = " + k.toPrecision(3) + " N/m";
		return k;
	}
	else{
		x = Math.sqrt(springEnergy*2/k);
		x = "\u0394x = " + x.toPrecision(3) + " m";
		return x;
	}
}
function energyGrav(){
	var gravEnergy = getData('gravener');
	var mass = getData('enermass');
	var height = getData('height');
	if(isNaN(gravEnergy)){
		gravEnergy = mass*9.8*height;
		gravEnergy = "Ug = " + gravEnergy.toPrecision(3) + " J";
		return gravEnergy;
	}
	else if(isNaN(mass)){
		mass = gravEnergy/(9.8*height);
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else{
		height = gravEnergy/(9.8*mass);
		height = "\u0394y = " + height.toPrecision(3) + " m";
		return height;
	}
}
function chooseEnerFormula(){
	var which = document.getElementById("energyselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "energyone"){
		answer = kineticEnergy();
	}
	else if(formula === "energytwo"){
		answer = calcWork();	
	}
	else if(formula === "energythree"){
		answer = calcPower();
	}
	else if(formula === "energyfour"){
		answer = rotKinetic();
	}
	else if(formula === "energyfive"){
		answer = energySpring();
	}
	else{
		answer = energyGrav();
	}
	document.getElementById("eneranswer").innerHTML = answer;
}
function changeEnergyVisible(){
	var which = document.getElementById("energyselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["kinenerpara", "enermasspara", "enervelpara", "workpara", "workforcepara", "distancepara", "powerpara", "enerchangepara","enertimepara", "rotkinpara", "rotinertpara", "eneromegapara", "enersprpara", "enerkpara", "enerstretchpara", "gravenerpara", "heightpara"];
	var ids2 = ["kinener", "enermass", "enervel", "work", "workforce", "distance", "power", "enerchange", "enertime", "rotkin", "rotinert", "eneromega", "enerspr", "enerk", "enerstretch", "gravener", "height"];
	var settings1 = ["inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings3 = ["none", "none", "none", "none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings4 = ["none", "none", "none", "none", "none", "none", "none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none"];
	var settings5 = ["none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "inline", "inline", "inline", "none", "none"];
	var settings6 = ["none", "inline", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "none", "inline", "inline"];
	if(formula === "energyone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings1[i];
		}
	}
	else if(formula === "energytwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];	
		}
	}
	else if(formula === "energythree"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];	
		}
	}
	else if(formula === "energyfour"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings4[i];	
		}
	}
	else if(formula === "energyfive"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings5[i];	
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings6[i];	
		}
	}	
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";	
	}
	document.getElementById("eneranswer").innerHTML = "";
}
function torqueNet(){
	var netTorque = getData('tnet');
	var rotInert = getData('inertia');
	var alpha = getData('alp');	
	if(isNaN(netTorque)){
		netTorque = rotInert*alpha;
		netTorque = "\u03C4net = " + netTorque.toPrecision(3) + " kgm^2/s^2";	
		return netTorque;
	}
	else if(isNaN(rotInert)){
		rotInert = netTorque/alpha;
		rotInert = "I = " + rotInert.toPrecision(3) + " kgm^2";
		return rotInert;	
	}
	else{
		alpha = netTorque/rotInert;
		alpha = "\u03B1 = " + alpha.toPrecision(3) + " rads/s^2";
		return alpha;
	}
}
function torque(){
	var t = getData('t');
	var radius = getData('trad');
	var force = getData('tforce');
	if(isNaN(t)){
		t = radius*force;
		t = "\u03C4 = " + t.toPrecision(3) + " kgm^2/s^2";
		return t;	
	}
	else if(isNaN(radius)){
		radius = t/force;
		radius = "r = " + radius.toPrecision(3) + " m";
		return radius;	
	}
	else{
		force = t/radius;
		force = "F = " + force.toPrecision(3) + " N";	
		return force;
	}
}
function angMoment(){
	var angularMom = getData('angmom');
	var inertia = getData('inertia');
	var omega = getData('momomega');
	if(isNaN(angularMom)){
		angularMom = inertia*omega;
		angularMom = "L = " + angularMom.toPrecision(3) + " 	kgm^2/s";
		return angularMom;
	}
	else if(isNaN(inertia)){
		inertia = angularMom/omega;
		inertia = "I = " + inertia.toPrecision(3) + " kgm^2";
		return inertia;	
	}
	else{
		omega = angularMom/inertia;
		omega = "\u03C9 = " + omega.toPrecision(3) + " rad/s";
		return omega;	
	}
}	
function changeAngMom(){
	var deltaAng = getData('delangmom');
	var t = getData('t');
	var time = getData('angmomtime');
	if(isNaN(deltaAng)){
		deltaAng = t*time;
		deltaAng = "\u0394L = " + deltaAng.toPrecision(3) + " kgm^2/s";
		return deltaAng;	
	}
	else if(isNaN(t)){
		t = deltaAng/time;
		t = "\u03C4 = " + t.toPrecision(3) + "kgm^2/s^2";
		return t;
	}
	else{
		time = deltaAng/t;
		time = "t = " + time.toPrecision(3) + " s";	
		return time;
	}
}
function chooseTorFormula(){
	var which = document.getElementById("torselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "torone"){
		answer = torqueNet();
	}
	else if(formula === "tortwo"){
		answer = torque();	
	}
	else if(formula === "torthree"){
		answer = angMoment();
	}
	else{
		answer = changeAngMom();
	}
	document.getElementById("toranswer").innerHTML = answer;
}
function changeTorVisible(){	
	var which = document.getElementById("torselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["tnetpara", "inertiapara", "alppara", "tpara", "tradpara", "tforcepara", "angmompara", "momomegapara", "delangmompara", "angmomtimepara"];
	var ids2 = ["tnet", "inertia", "alp", "t", "trad", "tforce", "angmom", "momomega", "delangmom", "angmomtime"];
	var settings = ["inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none"];
	var settings3 = ["none", "inline", "none", "none", "none", "none", "inline", "inline", "none", "none"];
	var settings4 = ["none", "none", "none", "inline", "none", "none", "none", "none", "inline", "inline"];
	if(formula === "torone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];	
		}
	}
	else if(formula === "tortwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];	
		}
	}
	else if(formula === "torthree"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];	
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings4[i];	
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";	
	}
	document.getElementById("toranswer").innerHTML = "";
}
function simpleHarmonicOne(){
	var x = getData('pos');
	var amplitude = getData('amplitude');
	var frequency = getData('freq');
	var time = getData('simpletime');
	if(isNaN(x)){
		x = amplitude*Math.cos(2*Math.PI*frequency*time);
		x = "x = " + x.toPrecision(3) + " m";
		return x;
	}
	else if(isNaN(amplitude)){
		amplitude = x/Math.cos(2*Math.PI*frequency*time);
		amplitude = "A = " + amplitude.toPrecision(3) + " m";
		return amplitude;
		
	}
	else if(isNaN(frequency)){
		frequency = Math.acos(x/amplitude)/(2*Math.PI*time);
		frequency = "f = " + frequency.toPrecision(3) + " 1/s";
		return frequency;
	}
	else{
		time = Math.acos(x/amplitude)/(2*Math.PI*frequency);
		time = "t = " + time.toPrecision(3) + " s";
		return time;
	}

}
function simpleHarmonicTwo(){
	var sprPeriod = getData('sprperiod');
	var mass = getData('simplemass');
	var k = getData('simplek');
	if(isNaN(sprPeriod)){
		sprPeriod = 2*Math.PI*Math.sqrt(mass/k);
		sprPeriod = "Ts = " + sprPeriod.toPrecision(3) + " s";
		return sprPeriod;
	}
	else if(isNaN(mass)){
		mass = Math.pow(sprPeriod/(2*Math.PI), 2)*k;
		mass = "m = " + mass.toPrecision(3) + " kg";
		return mass;
	}
	else{
		k = mass/Math.pow(sprPeriod/(2*Math.PI), 2);
		k = "k = " + k.toPrecision(3) + " N/m";
		return k;
	}
}
function simpleHarmonicThree(){
	var pendulumPeriod = getData('penperiod');
	var length = getData('length');
	var gravity = getData('gravity');
	if(isNaN(pendulumPeriod)){
		pendulumPeriod = 2*Math.PI*Math.sqrt(length/gravity);
		pendulumPeriod = "Tp = " + pendulumPeriod.toPrecision(3) + " s";
		return pendulumPeriod;
	}
	else if(isNaN(length)){
		length = Math.pow(pendulumPeriod/(2*Math.PI), 2)*gravity;
		length = "l = " + length.toPrecision(3) + " m";
		return length;
	}
	else{
		gravity = length/Math.pow(pendulumPeriod/(2*Math.PI), 2);
		gravity = "g = " + gravity.toPrecision(3) + " m/s^2";
		return gravity;
	}
}
function waves(){
	var f = getData('freq');
	var wavelength = getData('wavelength');
	var v = getData('wavevelocity');
	if(isNaN(f)){
		f = v/wavelength;
		f = "f = " + f.toPrecision(3) + " 1/s";
		return f;
	}
	else if(isNaN(wavelength)){
		wavelength = v/f;
		wavelength = "\u03BB = " + wavelength.toPrecision(3) + " m";
		return wavelength;
	}
	else{
		v = wavelength*f;
		v = "v = " + v.toPrecision(3) + " m/s";
		return v;
	}
}
function chooseSimpleFormula(){
	var which = document.getElementById("simpleselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "simpleone"){
		answer = simpleHarmonicOne();
	}
	else if(formula === "simpletwo"){
		answer = simpleHarmonicTwo();	
	}
	else if(formula === "simplethree"){
		answer = simpleHarmonicThree();
	}
	else{
		answer = waves();
	}
	document.getElementById("simpleanswer").innerHTML = answer;
}
function changeSimpleVisible(){
	var which = document.getElementById("simpleselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["pospara", "amplitudepara", "freqpara", "simpletimepara", "sprperiodpara", "simplemasspara", "simplekpara", "penperiodpara", "lengthpara", "gravitypara", "wavelengthpara", "wavevelocitypara"];
	var ids2 = ["pos", "amplitude", "freq", "simpletime", "sprperiod", "simplemass", "simplek", "penperiod", "length", "gravity", "wavelength", "wavevelocity"];
	var settings = ["inline", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none"];
	var settings3 = ["none", "none", "none", "none", "none", "none", "none", "inline" ,"inline", "inline", "none", "none"];
	var settings4 = ["none", "none", "inline", "none", "none", "none", "none", "none" ,"none", "none", "inline", "inline"];
	if(formula === "simpleone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];
		}
	}
	else if(formula === "simpletwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];
		}
	}
	else if(formula === "simplethree"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings4[i];
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";
	}
	document.getElementById("simpleanswer").innerHTML = "";
}
function elecForce(){
	var forceElec = getData('elecforce');
	var chargeOne = getData('chargeone');
	var chargeTwo = getData('chargetwo');
	var r = getData('elecrad');
	if(isNaN(forceElec)){
		forceElec = 9e9*chargeOne*chargeTwo/Math.pow(r, 2);
		forceElec = "Fe = " + forceElec.toPrecision(3) + " N";
		return forceElec;
	}
	else if(isNaN(chargeOne)){
		chargeOne = forceElec*Math.pow(r, 2)/(chargeTwo*9e9);
		chargeOne = "q1 = " + chargeOne.toPrecision(3) + " C";
		return chargeOne;	
	}
	else if(isNaN(chargeTwo)){
		chargeTwo = forceElec*Math.pow(r, 2)/(chargeOne*9e9);
		chargeTwo = "q2 = " + chargeTwo.toPrecision(3) + " C";
		return chargeTwo;	
	}
	else{
		r = Math.sqrt(chargeOne*chargeTwo*9e9/forceElec);
		r = "r = " + r.toPrecision(3) + " m";
		return r;	
	}
}
function elecTwo(){
	var I = getData('amps');
	var deltaCharge = getData('deltacharge');
	var t = getData('electime');
	if(isNaN(I)){
		I = deltaCharge/t;
		I = "I = " + I.toPrecision(3) + " A";
		return I;	
	}
	else if(isNaN(deltaCharge)){
		deltaCharge = I*t;
		deltaCharge = "	\u0394q = " + deltaCharge.toPrecision(3) + " C";
		return deltaCharge;	
	}
	else{
		time = deltaCharge/I;
		time = "t = " + time.toPrecision(3) + " s";
		return time;	
	}
}
function elecThree(){
	var R = getData('resistance');
	var row = getData('row');
	var length = getData('eleclength');
	var area = getData('area');	
	if(isNaN(R)){
		R = row*length/area;
		R = "R = " + R.toPrecision(3) + " \u03A9";
		return R;
	}
	else if(isNaN(row)){
		row = R*area/length;
		row = "\u03C1 = " + row.toPrecision(3) + " \u03A9m";
		return row;	
	}
	else if(isNaN(length)){
		length = R*area/row;
		length = "l = " + length.toPrecision(3) + " m";
		return length;	
	}
	else{
		area = row*length/R;
		area = "A = " + area.toPrecision(3) + " m^2";
		return area;	
	}
}
function elecFour(){
	var I = getData('amps');
	var V = getData('voltage');
	var R = getData('resistance');
	if(isNaN(I)){
		I = V/R;
		I = "I = " + I.toPrecision(3) + " A";
		return I;	
	}
	else if(isNaN(V)){
		V = I*R;
		V = "\u0394V = " + V.toPrecision(3) + " V";
		return V;	
	}
	else {
		R = V/I;
		R = "R = " + R.toPrecision(3) + " \u03A9";
		return R;	
	}
}
function elecFive(){
	var P = getData('elecpower');
	var V = getData('voltage');
	var I = getData('amps');
	if(isNaN(P)){
		P = I*V;
		P = "P = " + P.toPrecision(3) + " W";
		return P;	
	}
	else if(isNaN(V)){
		V = P/I;
		V = "\u0394V = " + V.toPrecision(3) + " V";
		return V;	
	}
	else{
		I = P/V;
		I = "I = " + I.toPrecision(3) + " A";
		return I;	
	}
}
function chooseElecFormula(){
	var which = document.getElementById("elecselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var answer;
	if(formula === "elecone"){
		answer = elecForce();
	}
	else if(formula === "electwo"){
		answer = elecTwo();	
	}
	else if(formula === "electhree"){
		answer = elecThree();
	}
	else if(formula === "elecfour"){
		answer = elecFour();
	}
	else{
		answer = elecFive();
	}
	document.getElementById("elecanswer").innerHTML = answer;
}
function changeElecVisible(){
	var which = document.getElementById("elecselect");
	var formula = which.options[which.selectedIndex].value.toString();
	var ids = ["elecforcepara", "chargeonepara", "chargetwopara", "elecradpara", "ampspara", "deltachargepara", "electimepara", "resistancepara", "rowpara", "eleclengthpara", "areapara", "voltagepara", "elecpowerpara"];
	var ids2 = ["elecforce", "chargeone", "chargetwo", "elecrad", "amps", "deltacharge", "electime", "resistance", "row", "eleclength", "area", "voltage", "elecpower"];
	var settings = ["inline", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none", "none", "none", "none"];
	var settings2 = ["none", "none", "none", "none", "inline", "inline", "inline", "none", "none", "none", "none", "none", "none"];
	var settings3 = ["none", "none", "none", "none", "none", "none", "none", "inline", "inline", "inline", "inline", "none", "none"];
	var settings4 = ["none", "none", "none", "none", "inline", "none", "none", "inline", "none", "none", "none", "inline", "none"];
	var settings5 = ["none", "none", "none", "none", "inline", "none", "none", "none", "none", "none", "none", "inline", "inline"];
	if(formula === "elecone"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings[i];	
		}
	}
	else if(formula === "electwo"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings2[i];	
		}
	}
	else if(formula === "electhree"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings3[i];	
		}
	}
	else if(formula === "elecfour"){
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings4[i];	
		}
	}
	else{
		for(var i in ids){
			document.getElementById(ids[i]).style.display = settings5[i];	
		}
	}
	for(var i in ids2){
		document.getElementById(ids2[i]).value = "";	
	}
	document.getElementById("elecanswer").innerHTML = "";
}