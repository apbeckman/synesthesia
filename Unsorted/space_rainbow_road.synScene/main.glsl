

			//******** Common Code Begins ********

mat3 m = mat3( 0.00,  0.80,  0.60,
              -0.80,  0.36, -0.48,
              -0.60, -0.48,  0.64 );
float hash( float n )
{
    return fract(sin(n)*43758.5453);
}

vec2 hash22( vec2 p )
{
    //p = mod(p, 4.0); // tile
    p = vec2(dot(p,vec2(175.1,311.7)),
             dot(p,vec2(260.5,752.3)));
    return fract(sin(p+455.)*18.5453);
}


float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*57.0 + 113.0*p.z;

    float res = mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
    return res;
}

float sdTriangle( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 )
{
	vec2 e0 = p1 - p0;
	vec2 e1 = p2 - p1;
	vec2 e2 = p0 - p2;

	vec2 v0 = p - p0;
	vec2 v1 = p - p1;
	vec2 v2 = p - p2;

	vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );
	vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
	vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );
    
    float s = e0.x*e2.y - e0.y*e2.x;
    vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),
                       vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),
                       vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));

	return -sqrt(d.x)*sign(d.y);
}

vec3 saturate(vec3 x)
{
    return clamp(x, vec3(0.0), vec3(1.0));
}


float sdCross(vec2 p) {
    p *= 2.;
    vec2 p2 = vec2(p.y, p.x);
    p2.y = abs(p2.y);
    p2.y -=.0;
    p.y = abs(p.y);
    p.y-=.0;

	vec2 v1 = vec2(-.5, 0.);
	vec2 v2 = vec2(.5, 0.);
	vec2 v3 = vec2(0., 1.5);
    float sd1= sdTriangle(p, v1, v2, v3);
    float sd2 = sdTriangle(p2, v1, v2, v3);
    return min(sd1, sd2);

}

 
mat2 rotate2d(float _angle){
    return mat2(cos(_angle),-sin(_angle),
                sin(_angle),cos(_angle));
}

float sdStar(vec2 p) {
    float d1 = sdCross(p);
    float d2 = sdCross(rotate2d(3.14 / 4.) * p  );
    return min(d1, d2);
}


float opExtrusion( in vec3 p, in float h )
{
    float d = sdStar(p.xy);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0));
}

			//******** BuffA Code Begins ********

float c1() {
    return (1. + clamp(-.2, .2, sin(TIME / 5.) ));
}

float c2() {
    return 0.;
    return .1*clamp(.0, 1., clamp(-.2, .2, sin(TIME / 3. + .7) ));
}

float hash31(vec3 p) {
   	float h = dot(p,vec3(127.1,311.7, 21.));	
    return fract(sin(h)*43758.5453123);
}

float hash21( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}

vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d )
{
    return a + b*cos( 6.28318*(c*t+d) );
    
}
 
vec3 HUEtoRGB(in float hue)
{
    vec3 rgb = abs(hue * 6. - vec3(3, 2, 4)) * vec3(1, -1, -1) + vec3(-1, 2, 2);
    return clamp(rgb, 0., 1.);
}

float path(float t) {
    return 0.1;
}

float multiwave(float path) {
    path /= 3.;
    return 10.*(2. * sin(path / 20.) + 
            1./4. * sin(path / 4.) +
            1./4. * sin(path / 4.));
}
float dmultiwave(float path) {
    path /= 3.;
    return (2. * cos(path / 20.) + 
            1./4. * cos(path / 4.) +
            1./4. * cos(path / 4.));
}

float opS( in float d1, in float d2 )
{
    d1 *= -1.0;
    return (d1 > d2) ? d1 : d2;
}


vec3 archCOL(vec3 p ) {
    float rep = 52.;
    float size = 11.;
    float in_size = 7.8;
    float l = 25.5;
    p.y += multiwave(p.z);
    float id = floor(p.z / rep);
    p.z = mod(p.z, rep) - rep / 2.;
    p.z-=20.;
    float outer = length(p.xy) - size;
    float iner = length(p.xy) - size + .5;
    float d = opS(iner, outer);
    
    float dist = max(d, abs( p.z - l) );
    vec3 c =  HUEtoRGB(hash31(vec3(id)));
    c /= 4.;
    c+= .3;
    return c;
}
 
float mapStarField(vec3 p, out vec2 id) {

    float r = c1();
    vec3 rep = vec3(20.);
    float scale = 1.;
    p.yz = rotate2d(3.14 / -2.) * p.yz;
    float t = TIME * 150.;
    p /= scale;
    p.y += t;
    float d = p.z;
    vec3 pm = mod(p, rep) - rep / 2.;
    id.x = hash31(floor(p / rep));
 
    if (id.x * r < 1.1) {
        return rep.z / 2.;
    }
    pm.xz = rotate2d( 42.*TIME *id.x) * pm.xz;
    id.x = fract(id.x *112.131);
    id.y = 1.;
   // id.x = .2;
    return opExtrusion(pm  / 2., .05);
}

float mapArch(vec3 p) {
    float rep = 52.;
    float size = 11.;
    float in_size = 7.8;
    float l = -25.5;
    p.y += multiwave(p.z);
    p.z = mod(p.z, rep) - rep / 2.;
    p.z-=20.;
    float outer = length(p.xy) - size;
    float iner = length(p.xy) - size + .5;
    float d = opS(iner, outer);
    
    return max(d, abs( p.z - l) );
}

float mapFence(vec3 p) {
    float rep = 2.;
    float height = .4;// + (1.+ sin(p.z / 4.)) / 1.;
    float size = .1 ;
    
    p.y += multiwave(p.z);
    p.z = mod(p.z, rep) - rep / 2.;
    p.x = abs(p.x) - 8.;
    float dt = length(vec2(p.x, p.y -.76 - height + 4.0*size)) - size / 2.;
        
    vec3 pstar = vec3(p.z, p.y, p.x);
    float dstar = opExtrusion(pstar, .01);
    float df = abs(max(length( p.xz - size), -.1+abs(p.y +.8) - height));
    return min(dt, dstar);
    
}

float mapRoad(vec3 p, out vec2 id) {
 
    
    vec2 rep = vec2(1, 1.);
    
    vec2 mp = mod(p.xz, rep) - rep / 2.;
    
    float s = .35;
    id.x = fract( ((+ p.z/rep.y)*0.323 + (p.x / rep.x )) *.0379523 );
   
    id.y = hash21( floor(p.xz / rep.xy ) * 10. );
    if (id.y < 0.02)
     id.x = fract( (floor(+ p.z/rep.y) *0.323+ floor(p.x / rep.x )) *0.0379523 );
    float l = TIME* 2. * hash21( floor(p.xz / rep)  );
    
     p.y -= .3 + sin(l ) * .0;
     float h =  p.y + 2. + multiwave(p.z) + cos(sin(p.z *c2()*.02  + TIME * c2() * .5) + p.x / 3.);
   
    return max( max(max(mp.x -s , mp.y - s ) , abs(p.x) - 8. ), abs( h )  - .05);
}

vec2 map(vec3 p, out vec2 id) {
    vec2 id2;
    vec3 ps = p;
    p.x+=sin(p.z / 15.)*2. + cos(p.z / 20.) + sin(p.z / 200.) * 200.;
    vec2 res =  vec2(mapRoad(p, id), 1.);
    
    vec2 resFence = vec2(mapFence(p), 2. );
    if (resFence.x < res.x) {
        res = resFence;
    }
    vec2 resArch = vec2(mapArch(p), 3.);
    if (resArch.x < res.x) {
        res = resArch;
    }
    vec2 resStar = vec2(mapStarField(ps, id2), 4.);
    if (resStar.x < res.x) {
        res = resStar;
        id = id2;
    }
    return res;
}

vec2 castRay(vec3 ro, vec3 rd, out vec2 id) {
    vec2 res = vec2(1e10, 0.);
    vec3 pos = ro;

    const float NEAR = 0.;
    const float FAR = 220.;

    float t = NEAR;
    
    for (int i = 0; i < 300 && t < FAR; ++i) {
    
        pos = ro + rd * t;
        vec2 h = map(pos, id);
        if (abs(h.x) < 0.001 * t) {
            res = vec2(t, h.y);
            break ;
        }
        t += t < 10. ? h.x * .5: h.x;
    }
    
    return res;
}


vec3 normal (in vec3 p)
{
    vec2 e = vec2(.0001, .0);
    vec2 id;
    float d = map (p,id).x;
    vec3 n = vec3 (map (p + e.xyy,id).x - d,
                   map (p + e.yxy,id).x - d,
                   map (p + e.yyx,id).x - d);
    return normalize(n);
}

vec3 raymarch_arch(vec3 ro, vec3 rd) {
    float t = 0.;
    vec3 acc = vec3(0.);

    for (int i = 0; i < 40; ++i) {
        vec3 pos = ro + rd * t;
         float dLight = .03+21.5/(10.+ length(mod(pos.zzz, 52.) -20.))*0.5;
        acc += dLight * archCOL(pos) / 50.;
    }
    return vec3(acc);
}

float stars(vec2 uv) {

    float acc = 0.;
    float s = 0.1;
    float scale = 100.;
    for (float i = 0.; i < 3.; i = i + 1.) {
        for (float j = 0.; j < 3.; j = j + 1.) {
            vec2 off = vec2( (j - 2.) * s, (i - 2.) * s);  
            acc += ( smoothstep( .92, 1.0, pow( noise(vec3( uv.x * scale + off.x, uv.y * scale + off.y, 0.)),3. )));
        }
    }
    return acc / 20.;
}

vec3 render(vec3 ro, vec3 rd, vec2 uv) {

    vec2 id;
    vec2 res = castRay(ro, rd, id);
    vec3 pos = ro + rd * res.x;
    float mat = res.y;
    float depth = length(pos - ro);
    vec3 albedo = vec3(20.);
    vec3 nor = normal(pos);
    vec3 lightDir = normalize( vec3(-2., 1., -0.));
    float ndotl = max(.4, dot(nor, lightDir));
   // vec3 stars = vec3( smoothstep( .9, 1.0, pow( noise(vec3( uv.x * 100., rd.y * 100., 0.)),3. )));
        vec3 stars = vec3(stars(vec2(uv.x, rd.y)));
    float dLight = .03+21.5/(5.+ length(mod(pos.zzz, 52.) -20.))*0.5;
     vec3 dLightCol = archCOL(pos);
    if (mat == 0.) {
        return stars;
    }
    if (mat == 1. || mat == 4.) {
        albedo = 4.* HUEtoRGB(id.x);// * (texture(iChannel1, pos.xz * 11.).x * 3. + .1);
        if (id.y < 0.02)
        albedo *= 32.;
        if (mat == 4.)
        albedo *= 12.;
     }
    if (mat == 2.) {
        albedo = vec3(0.5, 0.5, 0.2) * 20.;
    }
    float f = 1. / (1. +pow(depth * .0015 , 2.));

    vec3 r = raymarch_arch(ro, rd);
    float m = clamp(0., 1., 1.7 - f);
    return vec3(mix(dLightCol*dLight*albedo *ndotl + albedo / 30., stars  * 4., .8));
}

void camera(out vec3 ro, out vec3 rd, vec2 uv) {
    vec2 mouse = (-RENDERSIZE.xy + 2.0*_mouse.xy)/RENDERSIZE.y;
    float speed = TIME * 15.;
    ro = vec3( + sin(speed / 200.) * 200., 3.0 + multiwave(speed), -speed);
    mouse = vec2(0.);
    mouse.y += 3.;
    vec3 lookAt = vec3(sin(speed / 200. - 0.5 / 200.) * 200. + mouse.x + uv.x, mouse.y + -uv.y + multiwave(speed)- dmultiwave(speed)*.2, 1.-speed - c1() * 0.);
    
    rd = normalize(ro - lookAt);
    
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Normalized pixel coordinates (from -1 to 1)
    vec2 uv = (-RENDERSIZE.xy + 2.0*fragCoord)/RENDERSIZE.y;
    vec3 ro, rd;
    camera(ro, rd, uv);
    vec3 col = pow( render(ro, rd, uv)*1.6, vec3(1.0));
    fragColor = vec4(col * 1.,1.0);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

//First bloom pass, mipmap tree thing
float weights = 0.0;

vec3 ColorFetchB(vec2 coord)
{
 	return texture(BuffA, coord).rgb;   
}

vec3 Grab1(vec2 coord, const float octave, const vec2 offset)
{
 	float scale = exp2(octave);
    
    coord += offset;
    coord *= scale;

   	if (coord.x < 0.0 || coord.x > 1.0 || coord.y < 0.0 || coord.y > 1.0)
    {
     	return vec3(0.0);   
    }
    
    vec3 color = ColorFetchB(coord);

    return color;
}

vec3 Grab4(vec2 coord, const float octave, const vec2 offset)
{
 	float scale = exp2(octave);
    
    coord += offset;
    coord *= scale;

   	if (coord.x < 0.0 || coord.x > 1.0 || coord.y < 0.0 || coord.y > 1.0)
    {
     	return vec3(0.0);   
    }
    
    vec3 color = vec3(0.0);
    
    const int oversampling = 4;
    
    for (int i = 0; i < oversampling; i++)
    {    	    
        for (int j = 0; j < oversampling; j++)
        {
			vec2 off = (vec2(i, j) / RENDERSIZE.xy + vec2(0.0) / RENDERSIZE.xy) * scale / float(oversampling);
            color += ColorFetchB(coord + off);
            

            weights += 1.0;
        }
    }
    
    color /= weights;
    
    return color;
}

vec3 Grab8(vec2 coord, const float octave, const vec2 offset)
{
 	float scale = exp2(octave);
    
    coord += offset;
    coord *= scale;

   	if (coord.x < 0.0 || coord.x > 1.0 || coord.y < 0.0 || coord.y > 1.0)
    {
     	return vec3(0.0);   
    }
    
    vec3 color = vec3(0.0);
    float weights = 0.0;
    
    const int oversampling = 8;
    
    for (int i = 0; i < oversampling; i++)
    {    	    
        for (int j = 0; j < oversampling; j++)
        {
			vec2 off = (vec2(i, j) / RENDERSIZE.xy + vec2(0.0) / RENDERSIZE.xy) * scale / float(oversampling);
            color += ColorFetchB(coord + off);
            

            weights += 1.0;
        }
    }
    
    color /= weights;
    
    return color;
}

vec3 Grab16(vec2 coord, const float octave, const vec2 offset)
{
 	float scale = exp2(octave);
    
    coord += offset;
    coord *= scale;

   	if (coord.x < 0.0 || coord.x > 1.0 || coord.y < 0.0 || coord.y > 1.0)
    {
     	return vec3(0.0);   
    }
    
    vec3 color = vec3(0.0);
    float weights = 0.0;
    
    const int oversampling = 16;
    
    for (int i = 0; i < oversampling; i++)
    {    	    
        for (int j = 0; j < oversampling; j++)
        {
			vec2 off = (vec2(i, j) / RENDERSIZE.xy + vec2(0.0) / RENDERSIZE.xy) * scale / float(oversampling);
            color += ColorFetchB(coord + off);
            

            weights += 1.0;
        }
    }
    
    color /= weights;
    
    return color;
}

vec2 CalcOffset(float octave)
{
    vec2 offset = vec2(0.0);
    
    vec2 padding = vec2(10.0) / RENDERSIZE.xy;
    
    offset.x = -min(1.0, floor(octave / 3.0)) * (0.25 + padding.x);
    
    offset.y = -(1.0 - (1.0 / exp2(octave))) - padding.y * octave;

	offset.y += min(1.0, floor(octave / 3.0)) * 0.35;
    
 	return offset;   
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    
    vec3 color = vec3(0.0);
    
    /*
    Create a mipmap tree thingy with padding to prevent leaking bloom
   	
	Since there's no mipmaps for the previous buffer and the reduction process has to be done in one pass,
    oversampling is required for a proper result
	*/

    color += Grab1(uv, 1.0, vec2(0.0,  0.0)   );
    color += Grab4(uv, 2.0, vec2(CalcOffset(1.0))   );
    color += Grab8(uv, 3.0, vec2(CalcOffset(2.0))   );
    color += Grab16(uv, 4.0, vec2(CalcOffset(3.0))   );
    color += Grab16(uv, 5.0, vec2(CalcOffset(4.0))   );
    color += Grab16(uv, 6.0, vec2(CalcOffset(5.0))   );
    color += Grab16(uv, 7.0, vec2(CalcOffset(6.0))   );
    color += Grab16(uv, 8.0, vec2(CalcOffset(7.0))   );


    fragColor = vec4(color, 1.0);
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

//Horizontal gaussian blur leveraging hardware filtering for fewer texture lookups.

vec3 ColorFetchC(vec2 coord)
{
 	return texture(BuffB, coord).rgb;   
}
/*
float weights[5];
float offsets[5];
*/

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
  /*  
    weights[0] = 0.19638062;
    weights[1] = 0.29675293;
    weights[2] = 0.09442139;
    weights[3] = 0.01037598;
    weights[4] = 0.00025940;
    
    offsets[0] = 0.00000000;
    offsets[1] = 1.41176471;
    offsets[2] = 3.29411765;
    offsets[3] = 5.17647059;
    offsets[4] = 7.05882353;
    */
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec3 color = vec3(0.0);
    float weightSum = 0.0;
    
    if (uv.x < 1.12)
    {
        color += ColorFetchC(uv) * weights[0];
        weightSum += weights[0];

        for(int i = 1; i < 5; i++)
        {
            vec2 offset = vec2(offsets[i]) / RENDERSIZE.xy;
            color += ColorFetchC(uv + offset * vec2(0.5, 0.0)) * weights[i];
            color += ColorFetchC(uv - offset * vec2(0.5, 0.0)) * weights[i];
            weightSum += weights[i] * 2.0;
        }

        color /= weightSum;
    }

    fragColor = vec4(color,1.0);
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

//Vertical gaussian blur leveraging hardware filtering for fewer texture lookups.

vec3 ColorFetchD(vec2 coord)
{
 	return texture(BuffC, coord).rgb;   
}
/*
float weights[5];
float offsets[5];
*/

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
  
    weights[0] = 0.19638062;
    weights[1] = 0.29675293;
    weights[2] = 0.09442139;
    weights[3] = 0.01037598;
    weights[4] = 0.00025940;
    
    offsets[0] = 0.00000000;
    offsets[1] = 1.41176471;
    offsets[2] = 3.29411765;
    offsets[3] = 5.17647059;
    offsets[4] = 7.05882353;
    
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec3 color = vec3(0.0);
    float weightSum = 0.0;
    
    if (uv.x < 0.52)
    {
        color += ColorFetchD(uv) * weights[0];
        weightSum += weights[0];

        for(int i = 1; i < 5; i++)
        {
            vec2 offset = vec2(offsets[i]) / RENDERSIZE.xy;
            color += ColorFetchD(uv + offset * vec2(0.0, 0.5)) * weights[i];
            color += ColorFetchD(uv - offset * vec2(0.0, 0.5)) * weights[i];
            weightSum += weights[i] * 2.0;
        }

        color /= weightSum;
    }

    fragColor = vec4(color,1.0);
	return fragColor; 
 } 



vec4 cubic(float x)
{
    float x2 = x * x;
    float x3 = x2 * x;
    vec4 w;
    w.x =   -x3 + 3.0*x2 - 3.0*x + 1.0;
    w.y =  3.0*x3 - 6.0*x2       + 4.0;
    w.z = -3.0*x3 + 3.0*x2 + 3.0*x + 1.0;
    w.w =  x3;
    return w / 6.0;
}

vec4 BicubicTexture(in sampler2D tex, in vec2 coord)
{
	vec2 resolution = RENDERSIZE.xy;

	coord *= resolution;

	float fx = fract(coord.x);
    float fy = fract(coord.y);
    coord.x -= fx;
    coord.y -= fy;

    fx -= 0.5;
    fy -= 0.5;

    vec4 xcubic = cubic(fx);
    vec4 ycubic = cubic(fy);

    vec4 c = vec4(coord.x - 0.5, coord.x + 1.5, coord.y - 0.5, coord.y + 1.5);
    vec4 s = vec4(xcubic.x + xcubic.y, xcubic.z + xcubic.w, ycubic.x + ycubic.y, ycubic.z + ycubic.w);
    vec4 offset = c + vec4(xcubic.y, xcubic.w, ycubic.y, ycubic.w) / s;

    vec4 sample0 = texture(tex, vec2(offset.x, offset.z) / resolution);
    vec4 sample1 = texture(tex, vec2(offset.y, offset.z) / resolution);
    vec4 sample2 = texture(tex, vec2(offset.x, offset.w) / resolution);
    vec4 sample3 = texture(tex, vec2(offset.y, offset.w) / resolution);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix( mix(sample3, sample2, sx), mix(sample1, sample0, sx), sy);
}

vec3 ColorFetch(vec2 coord)
{
 	return texture(BuffA, coord).rgb;   
}

vec3 BloomFetch(vec2 coord)
{
 	return BicubicTexture(BuffD, coord).rgb;   
}

vec3 Grab(vec2 coord, const float octave, const vec2 offset)
{
 	float scale = exp2(octave);
    
    coord /= scale;
    coord -= offset;

    return BloomFetch(coord);
}

vec2 CalcOffset(float octave)
{
    vec2 offset = vec2(0.0);
    
    vec2 padding = vec2(10.0) / RENDERSIZE.xy;
    
    offset.x = -min(1.0, floor(octave / 3.0)) * (0.25 + padding.x);
    
    offset.y = -(1.0 - (1.0 / exp2(octave))) - padding.y * octave;

	offset.y += min(1.0, floor(octave / 3.0)) * 0.35;
    
 	return offset;   
}

vec3 GetBloom(vec2 coord)
{
 	vec3 bloom = vec3(0.0);
    
    //Reconstruct bloom from multiple blurred images
    bloom += Grab(coord, 1.0, vec2(CalcOffset(0.0))) * 1.0;
    bloom += Grab(coord, 2.0, vec2(CalcOffset(1.0))) * 1.5;
	bloom += Grab(coord, 3.0, vec2(CalcOffset(2.0))) * 1.0;
    bloom += Grab(coord, 4.0, vec2(CalcOffset(3.0))) * 1.5;
    bloom += Grab(coord, 5.0, vec2(CalcOffset(4.0))) * 1.8;
    bloom += Grab(coord, 6.0, vec2(CalcOffset(5.0))) * 1.0;
    bloom += Grab(coord, 7.0, vec2(CalcOffset(6.0))) * 1.0;
    bloom += Grab(coord, 8.0, vec2(CalcOffset(7.0))) * 1.0;

	return bloom;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec3 color = ColorFetch(uv);
    
    
    color += GetBloom(uv) * 0.12;
    
    color *= 2.0;
    

    //Tonemapping and color grading
    color = pow(color, vec3(1.5));
    color = color / (1.0 + color);
    color = pow(color, vec3(1.0 / 1.5));

    
    color = mix(color, color * color * (3.0 - 2.0 * color), vec3(1.0));
    color = pow(color, vec3(1.3, 1.20, 1.0));    

	color = saturate(color * 1.01);
    
    color = pow(color, vec3(0.7 / 2.2));

    fragColor = vec4(color, 1.0);
   // fragColor = texture(BuffB, uv) * 10111.;

	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}