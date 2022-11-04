vec4 iMouse = vec4(720, 720, 0.0, 0.0); 


 //    Transparent Isoslices
	// ---------------------

	// Quasi diffuse-lit transparent isosurface slices. I had to get a little inventive -
	// both with the character shuffling and the methodology - to make this work.

	// I don't know if anyone's noticed this, but 280 characters is not a lot to work 
	// with. :) I mean, you could order a beer in under two tweets, but you couldn't chat 
	// up the bar maid... OK, some people on Shadertoy probably could. :D

float udRoundBox( vec3 p, vec3 b, float r )
{
  return length(max(abs(p)-b,0.0))-r;
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

// Cheap organic distance field - partitioned into horizontal slices.
//
float m(vec3 p)
{ 
    // Moving the scene itself forward - as opposed to the camera - to save on some lighting
    // calculations later on. IQ does it in one of his small examples. "TIME" - IQ's suggestion.
    // I didn't know that was an option... I need to read more. :)
    // p.xy = _rotate(p.xy, p.z*roto_amt*sin(TIME+p.z));
    p.z += script_time; 

    
    // Chopping the field into horizontal strips. It's quite a brutal approach, and it only works
    // because of the volumetric-like traversal. By the way, if you comment the line out, the 
    // isolines disappear, along with 22 characters, but it wasn't the look I was after.
    p.y += 1/24.0;
    p.y = mix(p.y, ceil(p.y*12.0)/12.0, mix(sqrt(syn_OnBeat)*beat_slice, 1.0, always_sliced));
    // p.xy += 1/size_voxels*0.5;
    p = mix(p, ceil(0.5+p*size_voxels)/size_voxels, 1.0*voxelate);

    // Based on Public_int_i's suggestion. Saves a bunch of characters. If it matches the above
    // speedwise, I'll be using it.
    // p.y -= mod(p.y, .08); 
   
    // Very cheap organic distortion. This line can be amalgamated below... I was a little suspicious
    // about compilation order, but Greg Rostami is OK with it, so I am too. :)
    // p /= _triWave(p.z+bumpIn_scr*10.0-bumpOut_scr*10.0, 100.0);
    // p.z += syn_BPMTwitcher*0.1;
	// p = mix(p, floor(p)+smoothstep(vec3(0.0), vec3(1.0), fract(p)), ridges);
    // p *= _triWave(p.z+bumpIn_scr*10.0-bumpOut_scr*10.0, 100.0);
    // Return the absolute field value to ensure that the ray doesn't backtrack.
    float dist = abs(dot(p = cos(p*.6 + sin(p.zxy*1.8-morph_time*0.05)+sin(p.z*0.1+morph_time*0.1)), p) - 1.1);
    // dist = mix(dist, sdTorus(p, vec2(10.0))*0.1, ridges);
    return dist;
}

vec4 renderMain(){
    vec4 c = vec4(0.0);
    vec2 u = _xy;

    // The loosely centered and normalized unit direction ray, and the origin - initialized
    // to zero. Characters can be saved by dropping "d" directly into the loop, but this way
    // feels cleaner to me.
    vec3 d = vec3(u/RENDERSIZE.y - .6, 1), o = d - d;
    d.xz = _rotate(d.xz, camera_direction.x*0.7);
    d.yz = _rotate(d.yz, camera_direction.y);

    // Initializing "c" to zero... in a less than satisfactory way.
    c -= c;
    
	// Loop variables. Declaring the loop iteration variable, "i," outside of the loop in a 
    // shader environment makes me a little nervous, but things seem to be compiling.
    float j, i = 0.;
    
    // In order to save characters, it's customary to dispense with the early break.
    // Of course, that means we're effectively calling the function for the entire loop
    // count. Therefore, I figured I may as well try to take advantage of that... in as
    // few characters as possible. :)

    // Compiler's nightmare loop. :)
    //
    // Sketchy way to write a 100-count loop. Thanks to Aeikick for reminding me of that.
    // It's tempting to stick the first line after the second semicolon, but I think it'll be
    // skipped on the first iteration.
    //for(int i=99; i-->0;)
    // Slower "float" based loop, but uses fewer charaters. "i" is declared with "j"... It was 
    // BigWings's suggestion to  Nesvi7, so if it doesn't work, blame him. D
    //for(; i++<1e2;) 
    for (int i = 0; i<50; i++) // Cleaner Web 2.0 suggestion, but Abje.
        
        // The distance function. Normally, the "/3." wouldn't be there, but things have 
        // been shuffled around to save a few characters. By the way, the return distance, "m(o),"
        // is absolute to ensure that the ray continues travelling forward.    
        j = m(o)/3.,
        
        // If the distance to the surface comes within a certain threshold, add some lighting
        // based on said distance. The "m(o + const)" is an old trick. It's kind of a meld between 
        // directional derivative lighting and occlusion... Something like that anyway. "j" is
        // represetative of the current ray-to-surface distance. "o.z" is a rough approximation
        // of the total distance, and is used for attenuation. Using "o.z" has the added benefit of
        // giving the lighting a slight parallax feel.
        //
        // On a side note, using an "if" statement cuts down on characters, but I wanted this to
        // be branchless to save a few cycles.
        //
        // IQ's suggestion: Faster, and saves some characters. I'd originally used a step function,
        // which is apparently slower than the following. I read somewhere that the ternary 
        // operator isn't guaranteed not to branch on nontrivial calculations, but chances are
        // minimal, or something to that effect. I've also read that it never branches with certain
        // shader languages? ... I'll leave the debate to the experts. :)
        // Fabrice suggested a couple of obvious changes to cut it down some more.
        c -= j<.01 ? (j - .01)*m(o + .1)/(.02 + o.z/2e2) : 0.,
        // Advance the ray in such a way that it encourages it to pass closely to the surface, but
        // continue on its way.
        o += max(.01, j)*d;
                           

    // Cheap color... I'd need to cut this down further to fit it in, but I prefer the black and 
    // white look anyway.
    c.xyz += (d/o.z)*(syn_MidPresence+syn_MidLevel*syn_MidPresence+syn_MidLevel)*0.5*mids_color;
    c.zyx += (.1 + o*.2*sin(o.z))*colorize; //Etc.
    c += _loadImage(texIn, sin(_uv*10.0+TIME)*0.1+o.xy*0.2*(0.8+0.2*sin(TIME+_uvc.y*PI+_uv.x))+vec2(TIME*0.05, 0.0))*watery_tex;
    c += texture(texIn2, 0.5+o.xy*0.1*(0.8+0.2*sin(TIME+_uvc.y*PI+_uv.x))+vec2(TIME*0.05, 0.0))*digital_tex;
    // c += _loadImage(texIn, o.xy*0.2);

    if (syn_MediaType > 0.5){
        c += texture(syn_UserImage, 0.5+o.xy*0.1*(0.8+0.2*sin(TIME+_uvc.y*PI+_uv.x))+vec2(TIME*0.05, 0.0));
    }

    c += pow(c, vec4(3.0))*syn_HighHits*50.0*_pulse(sin(o.z*1.0-syn_HighTime), 1.0, 0.1)*flashing;
    // c *= (0.25+syn_FadeInOut*0.75);
	return c; 
 } 
