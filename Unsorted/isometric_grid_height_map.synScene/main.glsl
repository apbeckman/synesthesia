vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


/*

    Isometric Grid Height Map
    -------------------------
    
    Rendering some overlapping hexagons in an isometric grid fashion 
    to produce a height map of extruded blocks. Isometric grid 
    renderings are certainly not new, but they're fun and simple to 
    code.
    
    I had a few versions of this lying around for ages, but wasn't 
    inspired to do anything with them until Bitless posted his really 
    cool "Cyberspace data warehouse" example. He put isometric height
    maps on the faces of isometric cubes, which the isometric design 
    crowd have probably done before, but it was new to me. For anyone 
    who hasn't seen it, the link is below.
    
    A lot of people are familiar with the isometric height map concept, 
    but for anyone who isn't, you render cubes (or an extruded version) 
    in back to front order in a diamond pattern, or corresponding 
    patterns to cover the amount of overlap. In this case, I'm rendering 
    seven cells in a hexagon pattern.
    
    Just to mix things up a little and break the visual monotony, I
    rendered columns of differing width and offset them a bit. Other 
    than that, there's not a lot to it.
    
    By the way, rendering in this way is fun, novel and fast, but there 
    are way better 3D methods that will achieve the same.  I've attached 
    a link to a much simpler version of this for anyone interested in 
    the process, or code golfing an isometric heightmap, maze, etc. :)
    
    
    
    Inspired by the following:
    
    // Very cool, and unique.
    Cyberspace data warehouse - bitless
    https://www.shadertoy.com/view/NlK3Wt
    
    // Here's a more sophisticated example. It requires more
    // work, but is based on a similar concept.
    Isometric City 2.5D - knarkowicz
    https://www.shadertoy.com/view/MljXzz
    
    // A less interesting, but easier to understand, verion:
    Isometric Height Map - Shane
    https://www.shadertoy.com/view/flcSzX
    
*/

// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.619, 57.583)))*43758.5453); }


// Using two hexagons to creat an extruded box of varying height.
// I hacked in some logic to make it happen, but there might be
// more efficient ways. Either way, it works, so it'll do.
float dist(vec2 p, float h, float sz){

    //return length(p);
    
     // The vertical height component of the hexagon. We're taking
    // a little off to get rid of artifacts.
    float szD = sz/.8660254 - .015;
    
    // Height factor. More height looks more interesting, but too much
    // sends the blocks out of the cell range and causes artifacts.
    float hFact = .125+(syn_Intensity*0.45+syn_Presence*0.275);
    
    // Moving the blocks down a bit to maximize the range. I've hack a 
    // number in, but you could probably calculate something.
    p.y -= .1;
    
    
    // I should be rendering three calculated polygons consisting of
    // the top and two sides, but that requires more effort and GPU
    // power, so I've merely merged a base hexagon and a hexagon on
    // top that moves in the Y direction according to the height,
    // which is good enough for this example.
    vec2 q = abs(p);
    float hxB = max(q.y*.8660254 + q.x*.5, q.x);
    
    q = abs(p - vec2(0, szD - h)*hFact);
    float hxH = max(q.y*.8660254 + q.x*.5, q.x);

    // If the block is less that the size of a hexagon, you need the
    // maximum overlay of the two. If it's above the height of the
    // base hexagon, fuse the two together... It took me while to 
    // figure that out, but it works.
    float d = max(hxB, hxH);
    if(h - szD<0.) d = min(hxB, hxH);
    return d;
    
}

// If you're rendering triangles, hexagons, isometric grids, etc, the number
// 1.732 (sqrt(3)) will come up all the time.
const vec2 s = vec2(1.7320508, 1);

// Hexagon grid.
vec4 getGrid(vec2 p){

    // Finding the local coordinates and ID of the nearest hexagon center.
    vec4 ip = floor(vec4(p/s, p/s + .5));
    vec4 q = p.xyxy - vec4(ip.xy + .5, ip.zw)*s.xyxy;
    return dot(q.xy, q.xy)<dot(q.zw, q.zw)? vec4(q.xy, ip.xy) : vec4(q.zw, ip.zw - .5);
   
}

float height(vec2 p){

    //float rnd = hash21(p);
    //return smoothstep(.5, .95, sin(rnd*6.2831 + TIME)*.5 + .5);
    float tm = mod(smoothTimeC*0.5, 6.2831);
    float h = dot(sin(p*.73/1.25 - cos(p.yx*.97/1.25 - tm)*2.3), vec2(.25)) + .5;
    p *= 2.;
    float h2 = dot(sin(p*.73/1.25 - cos(p.yx*.97/1.25 - tm*2.)*2.3), vec2(.25)) + .5;
    h = mix(h, h2, .333);
    return smoothstep(0., 1., h);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;



    // Aspect correct coordinates.
    vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;

    // Scale and smoothing factor.
    const float sc = 16.;
    float sf = sc/RENDERSIZE.y;
    
    
    // Scaling and translation.
    vec2 p = uv*sc + s*smoothTime*0.6;
    // Including rotation for a side wall perpective.
    //vec2 p = rot2(-3.14159/2.)*uv*vec2(-1, 1)*sc + s.yx*TIME/4.;
    
    
    // Scene field calculations.
    
    // Light direction. Shining down and to the left.
    vec2 ld = normalize(s);
    
  
   
    // Rendering in a diamond grouping. Top first, the two below that, and the bottom.
    //vec2[4] cntr = vec2[4](vec2(0, .5), vec2(-.5, 0), vec2(.5, 0), vec2(0, -.5));
    
    
    // Rendering in a seven cell grouping, and taking rendering order into account:
    // Top cell first, then the two below, the one below those, two more below
    // that and the one on the bottom.
    vec2[7] cntr = vec2[7](vec2(0, .5), vec2(-.5, 0), vec2(.5, 0), vec2(0, -.5),
                           vec2(-.5, -.5), vec2(.5, -.5), vec2(0, -1));    
    
    float obj[7], objSh[7], side[7], side2[7], top[7], hgt[7];
    vec4[7] p4I;

    
    // The block width. The maximum is s.y/2. (.8660254) and will produce a packed 
    // grid with no gaps. However, due to overlap, the height needs to be reduced.
    float sz = .7;
    vec2 szOffs = sz/s.yx;
    
    // Object shadow.
    float objShad = 1e5;
    
    for(int i = min(0, int(FRAMECOUNT)); i<7; i++){
   
        // Obtain the hexagon and ID for this position. I could probably 
        // streamline the process, but this will do.
        vec4 p4 = getGrid(p + cntr[i]*s);
        // Random 2D offset.
        vec2 offs = vec2(hash21(p4.zw + .1), hash21(p4.zw + .2)) - .5;
        // Random size.
        sz = .7*(hash21(p4.zw + .3)*.35 + .65);
        szOffs = sz/s.yx; 
        // Position based height. 
        float h = height(p4.zw*s);
        vec2 q = p4.xy - cntr[i]*s - offs*vec2(.2, .2); 
        
        // The extruded column object, which is just two overlayed 
        // hexagons in disguise.
        obj[i] = dist(q, h, sz) - sz;
        
        /*
        // Failed experiment with different tops.
        float bx = -1e5;
        if(hash21(p4.zw + .23)<.5){
           bx = max(abs(q.x) - sz, 
                    abs(q.y - (sz/.8660254 - .015 - h)*.75 - .15) - sz*.5);
        }
        */
        
        // Using the column object above and some simple CSG to construt 
        // the top and sides.
        float sd1 = dist(q - vec2(-1)*szOffs, h, sz) - sz;
        float sd2 = dist(q - vec2(1, -1)*szOffs, h, sz) - sz;
        side[i] = max(obj[i], sd1);
        side2[i] = max(obj[i], sd2);
        top[i] = max(obj[i], -min(sd1, sd2));
        
        // Failed experiment with different tops.
        //top[i] = max(top[i] - .05, bx);  
        
        // The shadow of the column object.
        objSh[i] = dist(q - ld*.5*sz, h, sz) - sz;
        
        // Saving the height, and hexagon cell information.
        hgt[i] = h;
        p4I[i] = p4;
        
        // The shadow is laid down first on the background, so can be
        // calculated here.
        objShad = min(objShad,  objSh[i]);
 
    }
     
  

    
    // Initiate the scene color to the background.
    vec3 col = vec3(.143, .0125, .252);
    
 
    // Rendering the object shadows onto the background. 
    col = mix(col, vec3(.05, .05, .1), (1. - smoothstep(0., sf*2., objShad))*.4);
     
   
    // Loop through all seven object cells and rendering the objects.
    for(int i = min(0, int(FRAMECOUNT)); i<7; i++){
        
        // Dark glow around the objects for more fake AO.
        col = mix(col, vec3(0), (1. - smoothstep(0., sf*8., obj[i]))*.35);
        
        // Unique position-based cell ID.
        vec2 id = p4I[i].zw;
        
        
        // Coloring -- Mildly inspired by the colors in Shadertoy's in-house
        // "Rock Tiles" texture.
        float rnd = hash21(id);
        float fn = dot(sin(id*s/3.5/1.5 - cos(id.yx*s/2.3/1.5)*2.), vec2(.25)) + .5;
        fn = smoothstep(.1, .9, fn);  
        vec3 oCol = .6 + .4*cos(6.2831*mix(fn, rnd, .65)/3. + vec3(01.25, 2, 2.75)/.25);  
        oCol = mix(oCol, oCol.xzy, fn*(hash21(id + .3)*.25 + .25));
        if(hash21(id + .27)<.33) oCol = oCol*.6;
        if(hash21(id + .37)<.66) oCol = min(oCol*1.4, 1.);
        /*
        // Load in the "Rock Tiles" texture for a comparison.
        vec3 tx = texture(iChannel0, id*s/sc/1.5).xyz; tx *= tx;
        tx = smoothstep(-.05, .7, tx);
        vec3 oCol = tx*1.5;
        */
       
        // Using the height for some extra shadowing.
        float shad = 1. - hgt[i]*.5;
        
       
        // The top of the extruded box.
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, top[i])); // Edge, or stroke.
        //col = mix(col, oCol*ao*2., 1. - smoothstep(0., sf, top[i] + .04)); 
        //col = mix(col, vec3(0), 1. - smoothstep(0., sf, top[i] + .1)); 
        col = mix(col, oCol*shad, 1. - smoothstep(0., sf, top[i] + .04)); // Main.
        // Render dark holes on random faces to break up the monotony.
        if(hash21(id + .34)<.5){
        //if(mod(p4I[i].z + p4I[i].w, 2.)<.5){ // Checkered option.
            col = mix(col, vec3(0), 1. - smoothstep(0., sf, top[i] + .22)); 
            col = mix(col, oCol*shad/3., 1. - smoothstep(0., sf, top[i] + .265)); 
        }         

        
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, side[i])); // Edge, or stroke.
        //col = mix(col, vec3(.25, .5, .75)*shad, 1. - smoothstep(0., sf, side[i] + .04));
        //col = mix(col, vec3(0), 1. - smoothstep(0., sf, side[i] + .1));      
        col = mix(col, oCol*vec3(.25, .5, .75)*shad, 1. - smoothstep(0., sf, side[i] + .04)); // Main.
        
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, side2[i])); // Edge, or strke.
        //col = mix(col, vec3(2, .75, .5)*shad, 1. - smoothstep(0., sf, side2[i] + .04)); // Edge, or stroke.
        //col = mix(col, vec3(0), 1. - smoothstep(0., sf, side2[i] + .1)); // Edge, or strke.       
        col = mix(col, oCol*vec3(2, .75, .5)*shad, 1. - smoothstep(0., sf, side2[i] + .04)); // Main.
    
   
    }
   
    // Subtle vignette.
    uv = fragCoord/RENDERSIZE.xy;
    col *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .0625);

 
    // Rough gamma correction, and screen presentation.
    fragColor = vec4(sqrt(max(col, 0.)), 1);
    
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}