

			//******** BuffA Code Begins ********


vec3 glow = vec3(0);
#define dmin(a, b) a.x < b.x ? a : b
#define PI acos(-1.)
#define tau (2.*PI)
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define TIME (TIME + 3.6 + 5.)
#define pal(a,b,c,d,e) (a + b*sin(c*d + e))

vec3 att = vec3(1);

float pModPolar(inout vec2 p, float repetitions) {
	float angle = 2.*PI/repetitions;
	float a = atan(p.y, p.x) + angle/2.;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - angle/2.;
	p = vec2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2.)) c = abs(c);
	return c;
}
#define pmod(p,x) mod(p,x) - 0.5*x
vec4 valueNoise(float t){
	return mix(texture(image30,vec2(floor(t)/256.)),texture(image30,vec2(floor(t) + 1.)/256.), smoothstep(0.,1.,fract(t)));
}
// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2. * s)) - s)));
}

float iii;
float sdRhombus(vec3 p, vec3 s){
	
    
    p = abs(p) - s;
    
    float d = max(p.z, max(p.x, p.y));
    
    
    d = max(d, dot(p.yx + s.yx*0.5, normalize(vec2(1.))));
    d = max(d, dot(p.yz + s.yz*0.5, normalize(vec2(1.))));
    //d = max(d - s.x*0., -dot(p.z,p.x));
    //d 
    
    
    return d;
}

// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
	float e = max(r - abs(a - b), 0.);
	return min(a, b) - e*e*0.25/r;
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
	vec2 u = max(vec2(r - a,r - b), vec2(0));
	return max(r, min (a, b)) - length(u);
}
// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <r>):
float fOpUnionChamfer(float a, float b, float r) {
	return min(min(a, b), (a - r + b)*sqrt(0.5));
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	return c;
}

// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float fOpUnionColumns(float a, float b, float r, float n) {
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.)/((n-1.)*2.+sqrt(2.));
		pR45(p);
		p.x -= sqrt(2.)/2.*r;
		p.x += columnradius*sqrt(2.);
		if (mod(n,2.) == 1.) {
			p.y += columnradius;
		}
		// At this point, we have turned 45 degrees and moved at a point on the
		// diagonal that we want to place the columns on.
		// Now, repeat the domain along this direction and place a circle.
		pMod1(p.y, columnradius*2.);
		float result = length(p) - columnradius;
		result = min(result, p.x);
		result = min(result, a);
		return min(result, b);
	} else {
		return min(a, b);
	}
}



float sdVerticalCapsule( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}

float sdBox( vec3 p, vec3 s )
{
    p = abs(p) - s;
    return max(p.x, max(p.y, p.z));
}

    #define modD vec3(1.5,2.5,0.9)
float sdThangiPong(vec3 p){
    float mmm = sin(TIME);

    mmm = sign(mmm)*pow(abs(mmm), 5.);
    
    p.y -= TIME + iii + sin(TIME + iii);
    p.y = pmod(p.y, modD.y);
    
    
    p.xz *= rot(mmm*PI); 
    p.xz = abs(p.xz);
    float d = sdRhombus(p , vec3(0.12));
    
    glow += 0.5/(0.08+ d*d*4.)*vec3(0.1,0.43,0.3)*att;
    
    p.y -= 0.4;
    p.xz *= rot(0.25*PI);
    
    float n = fOpUnionStairs(d, sdRhombus(p - vec3(0,0.0 + 0.0,0), vec3(0.1)), 0.2,4.);
    d = min(d, n);
    //fOpUnionStairs
    return d;
}
float gg;
float speed = 0.5;

vec2 map(vec3 p){
	vec2 d = vec2(10e7);
    vec3 f = p;
    
    vec3 n = p;
    vec3 q = p;
    
    iii = pModPolar(q.xz, 3.);
    
    
    vec3 u = q;
    q.x-= 2.2;
    
    float dd = -q.x ;

    
    vec3 i = q;
    
    
    i.y = pmod(i.y, modD.y);
    vec3 k = i;
    i = abs(i);
    
    i = abs(i) - vec3(0.1,1.,0.4);
    float dm = max(i.z,max(i.y,i.x)) ;
    
    
    
    //glow += 0.1/(0.04 + dd*dd*150.)*vec3(1.,0.2,0.1);    
    
    
    q = abs(q);
    
    
    q.z -= 1.8;
    
    //float ddC = length(q.xz) - 0.1;
    float ddC = max(abs(q.x),abs(q.z)) - 0.2;
    
    
    dd = fOpUnionStairs(ddC, dd, 0.7,5.);
    
    
    
    q = pmod(q, modD);
    
    
    u = abs(q);
    float dW = min( min(
        			max(u.z, u.y),
        			max(u.x, u.y)
                ),        			max(u.x, u.z)) - 0.04 ;
    
    
    dd = fOpUnionStairs(dd, dW, 0.4,4.);
    
    
    
    float dF = p.y + 0.4;
    

    
    d = dmin(d, vec2(dd, 1.));
    
    d = dmin(d, vec2(dm, 11.));
    
    d = dmin(d, vec2(ddC, 1.));
    
    if(d.x == dW){
    	d.y = 6.;
    }
    
    n.y -= speed*TIME;
    
    

    
    float ddB = length(n) - 0.2;
    
    
    float mmm = pow(abs(sin(f.y*0.5 + TIME*0.5 + cos(f.z + TIME*0.25) +  + cos(f.y + TIME*0.125) )), 5.);
    float mmb = pow(abs(sin(f.y + TIME + cos(f.z*0.4 + TIME*0.45) +  + cos(f.y + TIME*0.5) )), 2.);
    
    float a= 0. + mmm*0.6;
    
    //ddB = fOpUnionStairs(ddB,ddD,2.,0.1 );
    
    
    n.xz *= rot(mmm*4.);
    n.yz *= rot(mmb*1.);
    
	vec3 j = n;
    
    j = pmod(j, 0.1);
    
    
    j = abs(j) - 0.06;
    float ddD = max(j.x, max(j.y, j.z));    
    
    n = abs(n);
    n.yz *= rot(0.9);
    n.xz -= 0.1 + mmb*0.1;
    
 	float dRr = sdRhombus(n , -vec3(0.04));
 	float dRb = sdRhombus(n , -vec3(0,0.1,0));
    
    
    //ddB = mix(ddB,dRr,mmb*0.5 );
    ddB = fOpUnionStairs(ddB,dRr,0.2 + mmb*0.1, 5. + sin(TIME) );
    ddB = fOpUnionStairs(ddB,dRb,0.25 + sin(TIME)*0.2, 5. );
    
    ddB = mix(ddB,ddD,a );
    d = dmin(d, vec2(ddB, 10.));
    
    
    k.x += 1.2;
    
    //k.y = abs(k.y);    
    
    
    
 	float dT = sdThangiPong(k);
    
    
    float dPp = length(k.xz) - 0.1; // pp huehue
    //glow += 0.1/(0.04 + dPp*dPp*10.)*vec3(1.,0.2,0.1)*att;    

    
    k.x -= 0.5;
    k = abs(k);
    k.z -= 0.9;
       
    //float dPpb = length(k.xz) - 0.1; // pp huehue
    
    k = abs(k);
    float dPpb = max(k.x,k.z) - 0.06; // pp huehue
    //glow += 0.1/(0.01 + dPpb*dPpb*14.)*vec3(0.6,0.1,0.1)*att;    
    glow += 0.2/(0.06 + dPpb*dPpb*dPpb*44.)*vec3(0.8,0.3,0.1)*att;    

    
    
    d = dmin(d, vec2(dT, 4.));
    d = dmin(d, vec2(abs(dPp) + 0.01, 4.));
    d = dmin(d, vec2(abs(dPpb) + 0.01, 4.));
    
    
    


    
    vec3 c = vec3(1);
	d.x *= 0.6;
    return d;
}
float dith;

float side;
vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = map(p);

    
    if(d.x < 0.3)
        ro += rd*0.3;
    p = ro; t = 0.; hit = false;
    for(int i = 0; i < 140; i++){
    	d = map(p);
        d.x *= dith * side;
        
    	//glow += exp(-d.x*20.);
        if(d.x < 0.001){
        	hit = true;
            break;
        }
        
        t += d.x;
        p = ro + rd*t;
    }
    
    
    return d;
}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
    vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0),dir ));
	vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up*uv.y);
}

vec3 getNormal(vec3 p){
	vec2 t= vec2(0.0004,0);  
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}

#define mx (2.*_mouse.x/RENDERSIZE.x)
#define my (0.6*_mouse.y/RENDERSIZE.x)
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    //
	//uv *= 1. - dot(uv,uv)*0.2;
    uv.xy *= rot(sin((TIME*0.7 - 3.6))*0.3);
    uv *= 1. + dot(uv,uv)*1.5;
    
    
    vec3 col = vec3(0);

    dith = mix(0.76,1., texture(image30, 20.*uv*256.).x);
    vec3 ro = vec3(0);
    
    ro.y += TIME*speed;
    ro.y -= 0.57 - my;
    
    float n = pow(valueNoise(TIME*1.).x, 2.);
    
    
    ro.y += sin(n);
    float nb = valueNoise(TIME*1./2.).x;
    float zoom = 1.2 + n*1.;
    n *= 1.;
    ro.xz += vec2(sin( nb*6.14*1.5 + mx),cos( nb*6.14*1.5 + mx))*zoom;
    
    ro.y += 0.3;
    
    vec3 lookAt = vec3(0,ro.y + sin(TIME)*0.05,0.);
    lookAt.y += -0.5 + valueNoise(TIME*1./2.).x;
    vec3 rd = getRd(ro, lookAt, uv);
    //rd.yz *= rot(TIME);
    
    vec3 p; float t = 0.; bool hit;
    float tA; side = 1.;
    
    for(int i = 0; i < 2   + min(0, FRAMECOUNT) ; i ++){
    	vec2 d = march(ro, rd, p, t, hit);
    	vec3 n = getNormal(p)*side;
        
        vec3 ld = normalize(vec3(1));
        vec3 j = normalize(ld - n);
        
        float diff = max(dot(n, ld), 0.);
        float ss = pow(max(dot(n, j), 0.), 20.);
        float fres = pow(1. - max(dot(n, -rd), 0.), 5.);
        
        //col += fres*0.04*vec3(0,0.5,1);
        //col += diff*fres*0.03*vec3(0.8,0.2,0.7);
        
        col += ss*0.05*vec3(1,0.1,1)*att;
        tA = max(tA,t);
        if (d.y == 10.){
            /*
        	side *= -1.;
            att *= vec3(0.2,0.6,1.)*0.9;
            
            rd = refract(rd, n,0.5);
            ro = p - n*0.4;*/
			
        	col += fres*0.1*vec3(0,0.5,1)*diff*att;
            ro = p + n*0.5;
            att *= vec3(0.2,0.7,1.)*(0.9 );
            
            //rd = -refract(rd, n,0.4);
            rd = reflect(rd, n);
        } else if (d.y == 11.){
            /*
        	side *= -1.;
            att *= vec3(0.2,0.6,1.)*0.9;
            
            rd = refract(rd, n,0.5);
            ro = p - n*0.1;*/
            
        	col += fres*0.1*vec3(0,0.5,1)*diff*att;
            ro = p + n*0.5;
            att *= vec3(0.2,0.6,1.)*0.9;   
            rd = reflect(rd, n + sin(p*40.)*0.01);
        }else if (d.y == 20.){
            /*
        	side *= -1.;
            att *= vec3(0.2,0.6,1.)*0.9;
            
            rd = refract(rd, n,0.5);
            ro = p - n*0.1;*/
            
        	col += fres*0.5*vec3(0.5,0.5,0.5)*glow*0.05*att;
            //ro = p + n*0.5;
            //att *= vec3(0.2,0.6,1.)*0.9;
            //rd = reflect(rd, n);
            break;
        } else {
            //col += fres*0.5*vec3(0.5,0.7,0.8)*0.3*att;
            #define aa(j) clamp(map(p + n/j).x*j, 0.,1.)
            //float aaa = aa(0.5)*aa(0.1)*aa(0.7)*20.;
            float aaa = aa(0.9)*aa(0.6)*10.;
            //aaa = 0.3;
            col += fres*0.5*vec3(0.5,0.7,0.8)*1.*att*aaa;
            //col += spec*0.5*vec3(0.5,0.7,0.4)*0.3*att;
        	break;
        }
    }
        
    col += glow*0.006;
    
    
    //col = mix(col, vec3(0.4,0.4,0.7)*0.3*att, pow(smoothstep(0.,1.,tA*0.143), 1.6));

    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 85" by jeyko. https://shadertoy.com/view/WdfczH
// 2020-03-14 09:43:07

// Fork of "Day 84" by jeyko. https://shadertoy.com/view/Wssczn
// 2020-03-13 10:52:05

// radial blur and chromatic abberation in this buffer
// thx iq for pallette and hg-sdf for polarMod


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 20.;
    float scale = 0.00 + pow(length(uv - 0.5),4.)*0.7;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*4.1;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor = radial*3.; 
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.6);
    //1fragColor *= 18.;
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*1.1));
    fragColor *= 1. - dot(uvn,uvn)*0.6;
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}