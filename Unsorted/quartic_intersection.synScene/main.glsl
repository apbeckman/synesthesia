

#define SHANE 1
#define MATTZ 2
#define WYATT 0
#define METHOD SHANE
mat3 M;
mat3 rot (vec3 u) {
		vec3 s = sin(u), c = cos(u);
		mat3 x = mat3(1,0,0, 		0,c.x,s.x, 		0,-s.x,c.x);
		mat3 y = mat3(c.y,0,s.y, 	0,1,0, 			-s.y,0,c.y);
		mat3 z = mat3(s.z,c.z,0,	-c.z,s.z,0,		0,0,1);
		return x*y*z;}

#if METHOD == WYATT
	//https://www.shadertoy.com/view/4dlcDN
	float cbrt (float a) {
        return sign(a)*exp(log(abs(a))/3.);
    }
    vec2 cbrt (vec2 a) {
        return vec2(cbrt(a.x),cbrt(a.y));
    }
    vec2 cbrti (vec2 v) {
        float a = length(v);
        a = exp(log(a)/3.);
        float theta = atan(v.y/v.x)/3.;
        return a*vec2(cos(theta), sin(theta));
    }
    float cubic (float b, float c, float d) {
        float p = -b/3.;
        float q = p*p*p + b*c/6. - 0.5*d;
        float r = c/3.;
        float s = r-p*p;
        float u = q*q + s*s*s;
        if (u < 0.) {
            vec2 t = vec2(0., sqrt(-u));
            vec2 Q = vec2(q, 0.);
            return (cbrti(Q + t) + cbrti(Q - t)).x + p; 
        } else {
            float t = sqrt(u);
            return cbrt(q + t) + cbrt(q - t) + p; 
        }
    }
	bvec4 solve_quartic (in vec4 coe, out vec4 roots) {
        //https://en.wikipedia.org/wiki/Quartic_function#Ferrari.27s_solution
        float b = coe.x;float c = coe.y; float d = coe.z; float e = coe.w;
        float p = c - 3.*b*b/8.;
        float q = 0.125*b*b*b - 0.5*b*c + d;
        float r = (-3.*b*b*b*b + 256.*e - 64.*b*d + 16.*b*b*c)/256.;
        float m = cubic(p, 0.25*p*p - r, -0.125*q*q);
        if (m<=0.) return bvec4(false);
        float A = -0.25*b;
        float B = 0.5*sqrt(2.*m);
        float C = -2.*p - 2.*m;
        float D = -sqrt(2.)*q/sqrt(m);
        float x = 10.;
        roots = vec4(
            A + B + 0.5*sqrt(C + D),
            A + B - 0.5*sqrt(C + D),
            A - B + 0.5*sqrt(C - D),
            A - B - 0.5*sqrt(C - D)
        );
        return bvec4(vec2(C+D>=0.),vec2(C-D>=0.));
    }

#elif METHOD == MATTZ
	//https://www.shadertoy.com/view/4dVcR1
	bvec4 solve_quartic(in vec4 coeffs,out vec4 roots) {
        //Mattz version of quartic solver
		//https://www.shadertoy.com/view/XdKyRR
        
    float p = coeffs[0];
    float q = coeffs[1]; 
    float r = coeffs[2];
    float s = coeffs[3];
    float i = -q;
    float j = p*r - 4.*s;
    float k = 4.*q*s - r*r - p*p*s;
    float a = (3.*j - i*i) / 3.;
    float b = (2.*i*i*i - 9.*i*j + 27.*k) / 27.;
    float delta1 = b*b / 4.;
    float delta2 = a*a*a / 27.;
    float delta = delta1 + delta2;
    
    float z1;
    
    if (delta >= 0.) {
        vec2 AB = -0.5*b + vec2(1,-1) * sqrt(max(delta, 0.));
        AB = sign(AB) * pow(abs(AB), vec2(1.0/3.0));
        z1 = AB.x + AB.y;
    } else {
        float phi = acos( -sign(b) * sqrt(delta1/-delta2) );
        z1 = 2. * sqrt(-a/3.) * cos( phi / 3.);
    }
    
    // shift back from normal form to root of resolvent cubic
    z1 -= i/3.;
    
    ////////////////////////////////////////////////////////////
	// now form quartic roots from resolvent cubic root

    float R2 = p*p/4. - q + z1; 
        
    bool R_ok = (R2 >= 0.);

    float R = sqrt(max(R2, 0.));
    
    float foo, bar;
    
    if (R == 0.) { 
        float z124s = z1*z1 - 4.*s;
        R_ok = R_ok && (z124s >= 0.);
        foo = 3.*p*p / 4. - 2.*q;
        bar = 2.*sqrt(max(z124s, 0.));
    } else {
        foo = 3.*p*p / 4. - R2 - 2.*q;
        bar = (4.*p*q - 8.*r - p*p*p) / (4.*R);
    }
    
    bool D_ok = R_ok && (foo + bar >= 0.);
    bool E_ok = R_ok && (foo - bar >= 0.);
    
    float D = sqrt(max(foo + bar, 0.));
    float E = sqrt(max(foo - bar, 0.));
    
    roots = vec4(-p/4.) + 0.5 * vec4(R+D, R-D, -(R-E), -(R+E));
    return bvec4(D_ok, D_ok, E_ok, E_ok);

}
#elif METHOD == SHANE
	//https://www.shadertoy.com/view/XsGyDh
	int solve_quadric(vec2 coeffs, inout vec2 roots){
        float p = coeffs.y / 2.;
        float D = p*p - coeffs.x;
        if (D <= 0.) return 0;
        else {
            roots = vec2(-1, 1)*sqrt(D) - p;
            return 2;
        }
    }
    int solve_cubic(vec3 coeffs, inout vec3 r){
        float a = coeffs[2];
        float b = coeffs[1];
        float c = coeffs[0];
        float p = b - a*a/3.;
        float q = a * (2.*a*a - 9.*b)/27. + c;
        float p3 = p*p*p;
        float d = q*q + 4.*p3/27.;
        float offset = -a/3.;
        if(d >= 0.0) { 
            vec2 uv = (vec2(1, -1)*sqrt(d) - q)/2.;
            uv = uv = sign(uv)*pow(abs(uv), vec2(1./3.));
            r[0] = offset + uv.x + uv.y;	
            float f = ((r[0] + a)*r[0] + b)*r[0] + c;
            float f1 = (3.*r[0] + 2. * a)*r[0] + b;
            r[0] -= f/f1;
            return 1;
        }
        float u = sqrt(-p/3.);
        float v = acos(-sqrt(-27./p3)*q/2.)/3.;
        float m = cos(v), n = sin(v)*1.732050808;
        float f,f1;
        r[0] = offset + u * (m + m);
        f = ((r[0] + a)*r[0] + b)*r[0] + c;
        f1 = (3.*r[0] + 2. * a)*r[0] + b;
        r[0] -= f / f1;
        r[1] = offset - u * (n + m);
        f = ((r[1] + a)*r[1] + b) * r[1] + c;
        f1=(3.*r[1] + 2. * a)*r[1] + b;
        r[1] -= f / f1;
        r[2] = offset + u * (n - m);
        f = ((r[2] + a)*r[2] + b)*r[2] + c;
        f1 = (3.*r[2] + 2. * a)*r[2] + b;
        r[2] -= f / f1;
        return 3;
    }
    bvec4 solve_quartic(vec4 coeffs, inout vec4 s){
        bvec4 broots;
        float a = coeffs[0];
        float b = coeffs[1];
        float c = coeffs[2];
        float d = coeffs[3];
        float sq_a = a * a;
        float p = - 3./8. * sq_a + b;
        float q = 1./8. * sq_a * a - 1./2. * a * b + c;
        float r = - 3./256.*sq_a*sq_a + 1./16.*sq_a*b - 1./4.*a*c + d;
        int num;
        vec3 cubic_coeffs;
        cubic_coeffs[0] = 1.0/2. * r * p - 1.0/8. * q * q;
        cubic_coeffs[1] = - r;
        cubic_coeffs[2] = - 1.0/2. * p;
        solve_cubic(cubic_coeffs, s.xyz);
        float z = s[0];
        float u = z * z - r;
        float v = 2. * z - p;
        if(u > 0.) u = sqrt(abs(u));
        else return bvec4(false);
        if(v > 0.) v = sqrt(abs(v));
        else return bvec4(false);
        vec2 quad_coeffs;
        quad_coeffs[0] = z - u;
        quad_coeffs[1] = q < 0. ? -v : v;
        num = solve_quadric(quad_coeffs, s.xy);
        if (num == 0) broots.xy = bvec2(false);
        if (num == 2) broots.xy = bvec2(true);
        quad_coeffs[0] = z + u;
        quad_coeffs[1] = q < 0. ? v : -v;
        vec2 tmp = vec2(1e8);
        int old_num = num;
        num = solve_quadric(quad_coeffs, s.zw);
        if (num == 0) broots.zw = bvec2(false);
        if (num == 2) broots.zw = bvec2(true);
        s -= a/4.;
        return broots;
    }
#endif

float absmin(float a, float b) {
	if (b>0.) return min(a,b);
    return a;
}
float intersect (vec4 coes) {
    vec4 roots;
    bvec4 br = solve_quartic(coes, roots);
	float i = 1e4;
    if (br.x) i = absmin(i,roots.x);
    if (br.y) i = absmin(i,roots.y);
    if (br.z) i = absmin(i,roots.z);
    if (br.w) i = absmin(i,roots.w);
    return i;
}
float intersect (vec2 coes) {
	float i=1e4;
    float det = coes.x*coes.x-4.*coes.y;
    if (det < 0.) return i;
    det =sqrt(det);
    i = absmin(i,0.5*(-coes.x+det));
    i = absmin(i,0.5*(-coes.x-det));
    return i;
}
vec4 torus (vec3 p, vec3 d, vec3 c, vec3 n, vec2 r) {
	float dn = dot(d,n);
    float wn = dot(p-c,n);
    vec3 s = p-c-wn*n;
    vec3 q = d - dn*n;
    float qq = dot(q,q);
    float sq = dot(s,q);
    float ss = dot(s,s);
    float A = (dn*dn+qq)*0.5/r.x;
    float B = (wn*dn+sq)/r.x;
    float C = (r.x*r.x-r.y*r.y+wn*wn+ss)*0.5/r.x;
    return vec4(2.*A*B,B*B+2.*A*C-qq,2.*C*B-2.*sq,C*C-ss)/(A*A);
    
}
vec4 cube (vec3 p, vec3 d, vec3 c, mat3 n, float r) {
	vec3 a = vec3(dot(p-c,n[0]),dot(p-c,n[1]),dot(p-c,n[2]));
    vec3 b = vec3(dot(d,n[0]),dot(d,n[1]),dot(d,n[2]));
    return vec4(
        4.*dot(a*b,b*b),
        6.*dot(a*a,b*b),
        4.*dot(a*a,a*b),
           dot(a*a,a*a)-r*r*r*r
    )/dot(b*b,b*b);
}
vec2 ellipse (vec3 p, vec3 d, vec3 a, vec3 b, float r) {
	a = p-a;b = p-b;
    float 
        rr = r*r,
        ad = dot(a,d),
        bd = dot(b,d),
        aa = dot(a,a),
        bb = dot(b,b);
    return vec2(
    	ad*aa-ad*bb+bd*bb-bd*aa-rr*(ad+bd),
        -aa*bb+0.25*(aa*aa+bb*bb+rr*rr)+0.5*(aa*bb-rr*(aa+bb))
    )/(ad*ad+bd*bd-rr-2.*ad*bd);
}
vec2 sphere (vec3 p, vec3 d, vec3 c, float r) {
	c = p-c;
    return vec2(2.*dot(c,d),dot(c,c)-r*r);
}
vec2 parabola (vec3 p, vec3 d, vec3 a, vec3 b) {
	vec3 n = normalize(b-a);
    a = p-a;
    b = p-b;
    float
        dn = dot(d,n),
        an = dot(a,n),
        bd = dot(b,d),
        bb = dot(b,b);
   	return vec2(
    	2.*(an*dn-bd),
        an*an-bb
    )/(dn*dn-1.);
}
float scene (vec3 p, vec3 d) {
    float i = 1e3;
    i = min(i,intersect(torus (p,d,vec3(0),vec3(0,0,1),vec2(0.15,0.02))));
    i = min(i,intersect(ellipse(p,d,vec3(-.1,-.08,0.04),vec3(0.1,-.08,-.04),.23)));
    i = min(i,intersect(cube (p,d,vec3(0),mat3(1,0,0,0,1,0,0,0,1),0.5)));
    i = min(i,intersect(cube (p,d,vec3(.1),mat3(1,0,0,0,1,0,0,0,1),0.08)));
  	i = min(i,intersect(sphere(p,d,vec3(0,.1,0),0.03)));
    i = min(i,intersect(parabola(p,d,vec3(-0.02,0,0),vec3(0))));
    return i;
    
}
vec3 nor (vec3 p, vec3 d, float i) {
	float e = 1e-4;
    vec3 t = normalize(cross(d,vec3(1,0,0)));
    vec3 a = p+e*t+d*scene(p+e*t,d);
    t = normalize(cross(d,t));
    vec3 b = p+e*t+d*scene(p+e*t,d);
    p = p+d*i;
	return normalize(cross(b-p,a-p));
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord/RENDERSIZE.xy*2.-1.)*RENDERSIZE.xy/RENDERSIZE.yy;
	M = rot(0.1243*smoothTime*vec3(0,vec2(1)+0.3*(_mouse.xy/RENDERSIZE.xy*2.-1.)));
    vec3 p = M*vec3(0.01*uv,-.3);
    vec3 d = M*normalize(vec3(uv,2));
    float i;
    vec3 n;
    vec3 col = vec3(0);
    for (int o = 0; o < 8; o++) {
    	i = scene(p,d);
        n = nor(p,d,i);
        p += i*d;
       	d = reflect(d,n);
        p += 0.001*d;
        col += (0.3+0.7*n)/float(2*o/3+2);
    }
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}