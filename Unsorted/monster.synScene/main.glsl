



float TK = 1.;

//float PI = 3.1415926535;



vec2 rot(vec2 p,float r){

	mat2 m = mat2(cos(r),sin(r),-sin(r),cos(r));

	return m*p;

}



vec2 pmod(vec2 p,float n){

	float np = 2.0*PI/n;

	float r = atan(p.x,p.y)-0.5*np;

	r = mod(r,np)-0.5*np;

	return length(p)*vec2(cos(r),sin(r));

}



float cube(vec3 p,vec3 s){

	vec3 q = abs(p);

	vec3 m = max(s-q,0.0);

	return length(max(q-s,0.0))-min(min(m.x,m.y),m.z);

}



float dist(vec3 p){
	p.xy = _rotate(p.xy,Spin*PI);
	p.xy = rot(p.xy,Twist*p.z);

	p.z -= TK*bass_time;


	p.xy = pmod(p.xy*Complexity,6.0);

	float k = 0.7;

	float zid = floor(p.z*k);

	p = mod(p,k)-0.5*k;

	for(int i = 0;i<6;i++){

		p = abs(p)-0.3;



		p.xy = rot(p.xy,1.0+zid+0.1*TK*spin_time*0.125);

		p.xz = rot(p.xz,1.0+4.7*zid+0.3*TK*spin_time*0.125);

	}

	return min(cube(p,vec3(0.3)),length(p)-0.4);

}





vec4 renderMainImage() {

	vec4 fragColor = vec4(0.0);

	vec2 fragCoord = _xy;



    vec2 uv = fragCoord/RENDERSIZE.xy;

	uv = 2.0*(uv-0.5);
	

	uv.y *= RENDERSIZE.y/RENDERSIZE.x;
	uv.xy += _uvc*FOV*PI;

	//			uv.xy += _uvc*FOV*PI;

	

	//uv = _rotate(uv, spin_time);

	uv += _rotate(uv, Rotate);

	vec3 ro = vec3(0.0,0.0,0.1);

	vec3 rd = normalize(vec3(uv+_uvc*PI*FOV,0.0)-ro);
	rd.xy = _rotate(rd.xy , spin_time);
	
	rd.xy*=1.0+Mirror*(rd.xy*PI*uv-1.0);

	rd.xz = _rotate(rd.xz, -PI*xyLook.x+PI*_uvc.x*Flip.x);
	rd.yz = _rotate(rd.yz, -PI*xyLook.y+PI*_uvc.y*Flip.y);

	float t  =2.0;

	float d = 0.0;

	float ac = 0.0;

	for(int i = 0;i<66;i++){

		d = dist(ro+rd*t)*0.12;

		d = max(0.0000,abs(d));

		t += d;

		if(d<0.001)ac += 0.1;//exp(-15.0*d);

	}

	vec3 col = vec3(0.0);

	col = vec3(0.21,0.27,0.537)*0.762*vec3(ac);//vec3(exp(-1.0*t));

	vec3 pn = ro+rd*t;

	float kn = 0.5;

	pn.z += -smoothTimeB*0.1;

	pn.z = mod(pn.z,kn)-0.5*kn;

	float em = clamp(0.01/pn.z,0.0,100.0);

//	col += 3.0*em*vec3(0.1,1.0,0.1)*pow(1.0+highhits*0.95, 2.);
	col += 3.0*em*vec3(01.1,.50,0.1);

	col = clamp(col,0.0,1.0);

	//col = 1.0-col;





    // Output to screen

    fragColor = vec4(col,0.);

	return fragColor; 

 } 





vec4 renderMain(){

	if(PASSINDEX == 0){

		return renderMainImage();

	}

}