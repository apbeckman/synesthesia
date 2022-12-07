

			//******** BuffA Code Begins ********

// Orientation
vec2 R;
vec2 hash(vec2 p) // Dave H
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);

}
vec3 rot (vec3 p, vec3 d) {
    float t = length(d);
    if (t==0.)return p;
    d = normalize(d);
    vec3 q = p-d*dot(d,p);
   	return p+(q)*(cos(t)-1.) + cross(d,q)*sin(t);
}
vec3 O (vec2 U) {return texture(BuffA,U/R).xyz;}
vec3 W (vec2 U) {return texture(BuffB,U/R).zyz;}
vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
    C.xyz = rot(O(U),W(U));
 	if (length(C.xyz)==0.) C.xyz = vec3(0,0,1);
 	vec3 ne = 0.125*(
        O(U+vec2(0, 1))+
        O(U+vec2(0,-1))+
        O(U+vec2( 1,0))+
        O(U+vec2(-1,0))+
        O(U+vec2(1, 1))+
        O(U+vec2(1,-1))+
        O(U+vec2( 1,1))+
        O(U+vec2(-1,1))+
        O(U+vec2(-1, 1))+
        O(U+vec2(-1,-1))+
        O(U+vec2( 1,-1))+
        O(U+vec2(-1,-1))
        
    );
 	C.xyz = mix(C.xyz,ne,.01);
 	C.xyz = normalize(C.xyz)*0.05;
 
 if (FRAMECOUNT <= 1) {
     C = vec4(0,0,1,0);
 	if (length(U-0.5*R)<2.) C.xyz = vec3(-1,0.001,0);
 }
 

	return C; 
 } 


			//******** BuffB Code Begins ********

//Angular momentum
//vec2 R;
/*vec2 hash(vec2 p) // Dave H
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);

}*/
vec3 F (vec3 A, vec3 B) {
    A = A-B;
	return A/dot(A,A)/length(A);
}
vec3 T (vec3 a, vec3 b, vec3 A, vec3 B) {
	vec3
        ap = a+A,
        an = a-A,
        bp = b+B,
        bn = b-B,
        Fp = F(ap,bp)-F(ap,bn),
        Fn = F(an,bn)-F(an,bp);
    return cross(A,Fp)-cross(A,Fn);
}
//vec3 O (vec2 U) {return texture(BuffA,U/R).xyz;}
//vec3 W (vec2 U) {return texture(BuffB,U/R).zyz;}
vec3 P (vec2 U) {
    if (fract(0.5*U.y)<0.5) U.x -= 0.5;
    U.x *= 1.732050808;
    vec2 h = hash(floor(U))*2.-1.;
    float t = smoothTime+h.x*1000.;
    return vec3(U+.5*h+0.5*vec2(sin(t),-cos(t)),0);
}
vec3 S (vec2 U) {
	vec3 
        t = vec3(0),
        o = O(U);
    for (int x = -2; x <= 2; x++)
    for (int y = -2; y <= 2; y++)
    {if (x==0&&y==0) continue;
    	t += T(P(U),P(U+2.481*vec2(x,y)),o,O(U+1.1*vec2(x,y)));
    }
	return W(U) + t;
}
vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
    C.xyz = S(U);
 	vec3 ne = 0.25*(
        W(U+vec2(0, 1))+
        W(U+vec2(0,-1))+
        W(U+vec2( 1,0))+
        W(U+vec2(-1,0))
    );
 	C.xyz = mix(C.xyz,ne,.5);
 	if (length(C.xyz)>.2)C.xyz=normalize(C.xyz)*0.2;
 	if (FRAMECOUNT <= 1) C = vec4(0,0,0,0);
	return C; 
 } 


#define H 1.73205080757
//vec2 R;
float intersect (vec2 coes) {
    float det = coes.x*coes.x-4.*coes.y;
    if (det < 0.) return 1e4;
    return 0.5*(-coes.x-sqrt(det));
}
float sphere (vec3 p, vec3 d, vec3 c, float r) {
	c = p-c;
    return intersect(vec2(2.*dot(c,d),dot(c,c)-r*r));
}
float ellipse (vec3 p, vec3 d, vec3 a, vec3 b, float r) {
	a = p-a;b = p-b;
    float 
        rr = r*r,
        ad = dot(a,d),
        bd = dot(b,d),
        aa = dot(a,a),
        bb = dot(b,b);
    return intersect(vec2(
    	ad*aa-ad*bb+bd*bb-bd*aa-rr*(ad+bd),
        -aa*bb+0.25*(aa*aa+bb*bb+rr*rr)+0.5*(aa*bb-rr*(aa+bb))
    )/(ad*ad+bd*bd-rr-2.*ad*bd));
}
vec3 norEllipse (vec3 p, vec3 a, vec3 b) {
    return normalize(normalize(p-a)+normalize(p-b));
}
float plane (vec3 p, vec3 d) {
	return  dot(-p,vec3(0,0,1))/dot(d,vec3(0,0,1));
}
vec3 On (vec2 U) {return normalize(texture(BuffA,U/R).xyz);}
//vec3 W (vec2 U) {return texture(iChannel1,U/R).zyz;}
float D (vec3 p, vec3 d, vec2 v, float e, inout vec4 color) {
    
    vec2 U = floor(p.xy+0.5+v);
    if (fract(0.5*U.y)<0.5) U.x += 0.5;
    vec3 o = On(U);
    vec2 h = hash(U)*2.-1.;
    float t = smoothTime+h.x*100.;
    vec3 q = vec3(U+0.5*h+0.3*vec2(sin(t),-cos(t)),sin(t)),
         a = q + o,
         b = q - o;
    float f = ellipse(p,d,a,b,2.1+basshits);
    if (f < e) {
    	e = f;
        p += d*e;
        vec3 n = norEllipse(p,a,b);
        float m = dot(n,o);
        color = dot(n,vec3(0,0,-1))*vec4(o.xzy*0.5+0.5,1);
    }
    return e;
}
vec4 X (vec3 p, vec3 d) {
	p += d*plane(p, d);
	vec4 color = vec4(0);
    float e = 1e3;
    for (int x = -1; x<=1; x++)
    for (int y = -1; y<=1; y++)
        e = D (p, d, vec2(x,y), e, color);
    return color;
}
vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
 	U = 2.*(U-0.5*R)/R.y;
 	vec3 p = vec3(_mouse.xy+19.*U*Zoom,-1),
  
    d = normalize(vec3(U*0.1,1));
 	C = X(p,d);
 	vec3 o = On(gl_FragCoord.xy);
 	if (_mouse.z <1.) C.xyz = (o.xzy*0.5+0.5);
	return C; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderMainImage();
	}
}