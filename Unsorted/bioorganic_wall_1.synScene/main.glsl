

/*

	Bioorganic Wall
	---------------

	Raymarching a textured XY plane. Basically, just an excuse to try out the new pebbled texture.

*/

//	Classic Perlin 3D Noise 
//	by Stefan Gustavson
//
vec4 permute(vec4 x){return mod(((x*34.0)+1.0)*x, 289.0);}
vec4 taylorInvSqrt(vec4 r){return 1.79284291400159 - 0.85373472095314 * r;}
vec3 fade(vec3 t) {return t*t*t*(t*(t*6.0-15.0)+10.0);}

float cnoise(vec3 P){
  vec3 Pi0 = floor(P); // Integer part for indexing
  vec3 Pi1 = Pi0 + vec3(1.0); // Integer part + 1
  Pi0 = mod(Pi0, 289.0);
  Pi1 = mod(Pi1, 289.0);
  vec3 Pf0 = fract(P); // Fractional part for interpolation
  vec3 Pf1 = Pf0 - vec3(1.0); // Fractional part - 1.0
  vec4 ix = vec4(Pi0.x, Pi1.x, Pi0.x, Pi1.x);
  vec4 iy = vec4(Pi0.yy, Pi1.yy);
  vec4 iz0 = Pi0.zzzz;
  vec4 iz1 = Pi1.zzzz;

  vec4 ixy = permute(permute(ix) + iy);
  vec4 ixy0 = permute(ixy + iz0);
  vec4 ixy1 = permute(ixy + iz1);

  vec4 gx0 = ixy0 / 7.0;
  vec4 gy0 = fract(floor(gx0) / 7.0) - 0.5;
  gx0 = fract(gx0);
  vec4 gz0 = vec4(0.5) - abs(gx0) - abs(gy0);
  vec4 sz0 = step(gz0, vec4(0.0));
  gx0 -= sz0 * (step(0.0, gx0) - 0.5);
  gy0 -= sz0 * (step(0.0, gy0) - 0.5);

  vec4 gx1 = ixy1 / 7.0;
  vec4 gy1 = fract(floor(gx1) / 7.0) - 0.5;
  gx1 = fract(gx1);
  vec4 gz1 = vec4(0.5) - abs(gx1) - abs(gy1);
  vec4 sz1 = step(gz1, vec4(0.0));
  gx1 -= sz1 * (step(0.0, gx1) - 0.5);
  gy1 -= sz1 * (step(0.0, gy1) - 0.5);

  vec3 g000 = vec3(gx0.x,gy0.x,gz0.x);
  vec3 g100 = vec3(gx0.y,gy0.y,gz0.y);
  vec3 g010 = vec3(gx0.z,gy0.z,gz0.z);
  vec3 g110 = vec3(gx0.w,gy0.w,gz0.w);
  vec3 g001 = vec3(gx1.x,gy1.x,gz1.x);
  vec3 g101 = vec3(gx1.y,gy1.y,gz1.y);
  vec3 g011 = vec3(gx1.z,gy1.z,gz1.z);
  vec3 g111 = vec3(gx1.w,gy1.w,gz1.w);

  vec4 norm0 = taylorInvSqrt(vec4(dot(g000, g000), dot(g010, g010), dot(g100, g100), dot(g110, g110)));
  g000 *= norm0.x;
  g010 *= norm0.y;
  g100 *= norm0.z;
  g110 *= norm0.w;
  vec4 norm1 = taylorInvSqrt(vec4(dot(g001, g001), dot(g011, g011), dot(g101, g101), dot(g111, g111)));
  g001 *= norm1.x;
  g011 *= norm1.y;
  g101 *= norm1.z;
  g111 *= norm1.w;

  float n000 = dot(g000, Pf0);
  float n100 = dot(g100, vec3(Pf1.x, Pf0.yz));
  float n010 = dot(g010, vec3(Pf0.x, Pf1.y, Pf0.z));
  float n110 = dot(g110, vec3(Pf1.xy, Pf0.z));
  float n001 = dot(g001, vec3(Pf0.xy, Pf1.z));
  float n101 = dot(g101, vec3(Pf1.x, Pf0.y, Pf1.z));
  float n011 = dot(g011, vec3(Pf0.x, Pf1.yz));
  float n111 = dot(g111, Pf1);

  vec3 fade_xyz = fade(Pf0);
  vec4 n_z = mix(vec4(n000, n100, n010, n110), vec4(n001, n101, n011, n111), fade_xyz.z);
  vec2 n_yz = mix(n_z.xy, n_z.zw, fade_xyz.y);
  float n_xyz = mix(n_yz.x, n_yz.y, fade_xyz.x); 
  return 2.2 * n_xyz;
}


// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ){
   
    n = max(n*n, 0.001);
    n /= (n.x + n.y + n.z );  
    
	return (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
    
}

// Raymarching a textured XY-plane, with a bit of distortion thrown in.
float map(vec3 p){
    float F = (cnoise(vec3(smoothTime*0.001)));
    // A bit of cheap, lame distortion for the heaving in and out effect.
    p.xy += sin(p.xy*7. + cos(p.yx*13. + smoothTimeC*0.25))*.01;
    
    // Back plane, placed at vec3(0., 0., 1.), with plane normal vec3(0., 0., -1).
    // Adding some height to the plane from the texture. Not much else to it.
    return 1. - p.z - texture(image12, p.xy+cnoise(F*0.1+cos(p.xyy+(smoothTimeC*0.001)))).x*.125;

    
    // Flattened tops.
    //float t = texture(image47, p.xy).x;
    //return 1. - p.z - smoothstep(0., .7, t)*.06 - t*t*.03;
    
}


// Tetrahedral normal, courtesy of IQ.
vec3 getNormal( in vec3 pos )
{
    vec2 e = vec2(0.001, -0.001);
    return normalize(
        e.xyy * map(pos + e.xyy) + 
        e.yyx * map(pos + e.yyx) + 
        e.yxy * map(pos + e.yxy) + 
        e.xxx * map(pos + e.xxx));
}

// Ambient occlusion, for that self shadowed look.
// Based on the original by IQ.
float calculateAO(vec3 p, vec3 n){

   const float AO_SAMPLES = 5.0;
   float r = 1.0, w = 1.0, d0;
    
   for (float i=1.0; i<=AO_SAMPLES; i++){
   
      d0 = i/AO_SAMPLES;
      r += w * (map(p + n * d0) - d0);
      w *= 0.5;
   }
   return clamp(r, 0.0, 1.0);
}

// Cool curve function, by Shadertoy user, Nimitz.
//
// It gives you a scalar curvature value for an object's signed distance function, which 
// is pretty handy for all kinds of things. Here, it's used to darken the crevices.
//
// From an intuitive sense, the function returns a weighted difference between a surface 
// value and some surrounding values - arranged in a simplex tetrahedral fashion for minimal
// calculations, I'm assuming. Almost common sense... almost. :)
//
// Original usage (I think?) - Cheap curvature: https://www.shadertoy.com/view/Xts3WM
// Other usage: Xyptonjtroz: https://www.shadertoy.com/view/4ts3z2
float curve(in vec3 p){

    const float eps = 0.02, amp = 7., ampInit = 0.5;

    vec2 e = vec2(-1., 1.)*eps; //0.05->3.5 - 0.04->5.5 - 0.03->10.->0.1->1.
    
    float t1 = map(p + e.yxx), t2 = map(p + e.xxy);
    float t3 = map(p + e.xyx), t4 = map(p + e.yyy);
    
    return clamp((t1 + t2 + t3 + t4 - 4.*map(p))*amp + ampInit, 0., 1.);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    
    // Unit directional ray.
    vec3 rd = normalize(vec3(fragCoord - RENDERSIZE.xy*.5, RENDERSIZE.y*1.5));
    
    rd.yz = _rotate(rd.yz,_uvc.y/PI);
    
    // Rotating the XY-plane back and forth, for a bit of variance.
    // Compact 2D rotation matrix, courtesy of Shadertoy user, "Fabrice Neyret."
    vec2 a = sin(vec2(1.5707963, 0) + sin(smoothTimeC*0.025)*.3);
    //rd.xy = mat2(a, -a.y, a.x)*rd.xy;
    // Ray origin. Moving in the X-direction to the right.
    vec3 ro = vec3(0.+sin(TIME*0.05), smoothTime*.035, 0.);
    
    
    // Light position, hovering around camera.
    vec3 lp = ro + vec3(cos(smoothTimeB*0.125)*.5, sin(smoothTimeB*0.125)*.5, 0.5);
    
    // Standard raymarching segment. Because of the straight forward setup, very few 
    // iterations are needed.
    float d, t=0.;
    for(int j=0;j<16;j++){
      
        d = map(ro + rd*t); // distance to the function.
        
        // The plane "is" the far plane, so no far plane break is needed.
        if(d<0.001) break; 
        
        t += d*.7; // Total distance from the camera to the surface.
    
    }
    
   
    // Surface postion, surface normal and light direction.
    vec3 sp = ro + rd*t;
    vec3 sn = getNormal(sp);
    vec3 ld = lp - sp;
    
    
    // Retrieving the texel at the surface postion. A tri-planar mapping method is used to
    // give a little extra dimension. The time component is responsible for the texture movement.
    float c = .95 - tex3D(image47, sp*5. - vec3(sp.x, sp.y+cnoise(sp.xyy+smoothTimeC*0.001), smoothTimeC*0.0125+sp.x+sp.y), sn).x;
    
    // Taking the original grey texel shade and colorizing it. Most of the folowing lines are
    // a mixture of theory and trial and error. There are so many ways to go about it.
    //
    vec3 orange = vec3(min(c*1.5, 1.), pow(c, 2.), pow(c, 8.)); // Cheap, orangey palette.
    
    vec3 oC = orange; // Initializing the object (bumpy wall) color.
    
    // Old trick to shift the colors around a bit. Most of the figures are trial and error.
    oC = mix(oC, oC.zxy, cos(rd.zxy*6.283 + sin(sp.yzx*6.283))*.25+.75);
    oC = mix(oC.yxz, oC, (sn)*.25+.25); // Using the normal to colorize.
    
    oC = mix(orange, oC, (sn)*.25+.75);
    oC *= oC*1.5;
    
    // Plain, old black and white. In some ways, I prefer it. Be sure to comment out the above, though.
    //vec3 oC = vec3(pow(c, 1.25)); 
    
    
    // Lighting.
    //
    float lDist = max(length(ld), 0.001); // Light distance.
    float atten = 1./(1. + lDist*.125); // Light attenuation.
    
    ld /= lDist; // Normalizing the light direction vector.
    
    float diff = max(dot(ld, sn), 0.); // Diffuse.
    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 128.); // Specular.
    float fre = clamp(dot(sn, rd) + .25, .0, 1.); // Fake fresnel, for the glow.

    
    // Shading.
    //
    // Note, there are no actual shadows. The camera is front on, so the following two 
    // functions are enough to give a shadowy appearance.
    float crv = curve(sp); // Curve value, to darken the crevices.
    float ao = calculateAO(sp, sn); // Ambient occlusion, for self shadowing.

    // Not all that necessary, but adds a bit of green to the crevice color to give a fake,
    // slimey appearance.
    //vec3 crvC = vec3(crv, crv*1.3, crv*.7)*.25 + crv*.75;
    //crvC *= crvC;
    
    // Combining the terms above to light up and colorize the texel.
    vec3 col = (oC*(diff + .5) + vec3(.5, .75, 1.)*spec*2.) + vec3(.3, .7, 1.)*pow(fre, 3.)*5.;
    // Applying the shades.
    //col *= (atten*crvC*ao);
    col *= (atten*ao);
    
    // I might be stating the obvious here, but if you ever want to see what effect each individual
    // component has on the overall scene, just put the variable in at the end, by itself.
    //col = vec3(ao); // col = vec3(crv); // etc.

    // Presenting to the screen.
	fragColor = vec4(sqrt(max(col, 0.)), 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}