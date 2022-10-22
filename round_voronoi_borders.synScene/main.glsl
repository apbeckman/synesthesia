//******** Common Code Begins ********

// The Shadertoy variable "TIME" isn't recognized in the common tab, so this is a
// hacky workaround. Take a look at the fBm function.
float gTime;


// Unsigned distance to the segment joining "a" and "b."
float distLine(vec2 a, vec2 b){

	b = a - b;
	float h = clamp(dot(a, b) / dot(b, b), 0., 1.);
    return length(a - b*h);
}


vec2 hash22(vec2 p, float repScale) {

    // Repetition.
    p = mod(p, repScale);

        // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked
    // amalgamation I put together, based on a couple of other random algorithms I've
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(113, 1)));

    return (fract(vec2(262144, 32768)*n) - .5)*2. - 1.;
}


// Standard 2x2 hash algorithm... I'll replace this with Dave Hoskins's more reliable
// one at some stage.
vec2 hash22G(vec2 p, float repScale) {

    // Repetition.
    p = mod(p, repScale);

        // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked
    // amalgamation I put together, based on a couple of other random algorithms I've
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(113, 1)));

    // Static.
    //return (fract(vec2(262144, 32768)*n) - .5)*2. - 1.;

    // Animated.
    p = fract(vec2(262144, 32768)*n);
    // Note the ".45," insted of ".5" that you'd expect to see. When edging, it can open
    // up the cells ever so slightly for a more even spread. In fact, lower numbers work
    // even better, but then the random movement would become too restricted. Zero would
    // give you square cells.
    return sin( p*6.2831853 + morph_time);
}

// Gradient noise: Ken Perlin came up with it, or a version of it. Either way, this is
// based on IQ's implementation. It's a pretty simple process: Break space into squares,
// attach random 2D vectors to each of the square's four vertices, then smoothly
// interpolate the space between them.
float gradN2D(in vec2 f, float repScale){

   f *= repScale;

    // Used as shorthand to write things like vec3(1, 0, 1) in the short form, e.yxy.
   const vec2 e = vec2(0, 1);

    // Set up the cubic grid.
    // Integer value - unique to each cube, and used as an ID to generate random vectors for the
    // cube vertiies. Note that vertices shared among the cubes have the save random vectors attributed
    // to them.
    vec2 p = floor(f);
    f -= p; // Fractional position within the cube.


    // Smoothing - for smooth interpolation. Use the last line see the difference.
    //vec2 w = f*f*f*(f*(f*6.-15.)+10.); // Quintic smoothing. Slower and more squarish, but derivatives are smooth too.
    vec2 w = f*f*(3. - 2.*f); // Cubic smoothing.
    //vec2 w = f*f*f; w = ( 7. + (w - 7. ) * f ) * w; // Super smooth, but less practical.
    //vec2 w = .5 - .5*cos(f*3.14159); // Cosinusoidal smoothing.
    //vec2 w = f; // No smoothing. Gives a blocky appearance.

    // Smoothly interpolating between the four verticies of the square. Due to the shared vertices between
    // grid squares, the result is blending of random values throughout the 2D space. By the way, the "dot"
    // operation makes most sense visually, but isn't the only metric possible.
    float c = mix(mix(dot(hash22G(p + e.xx, repScale), f - e.xx), dot(hash22G(p + e.yx, repScale), f - e.yx), w.x),
                  mix(dot(hash22G(p + e.xy, repScale), f - e.xy), dot(hash22G(p + e.yy, repScale), f - e.yy), w.x), w.y);

    // Taking the final result, and converting it to the zero to one range.
    return c*.5 + .5; // Range: [0, 1].
}

// Gradient noise fBm.
float fBm(in vec2 p, float repScale, float time){

    gTime = time;

    return gradN2D(p, repScale)*.57 + gradN2D(p, repScale*2.)*.28 + gradN2D(p, repScale*4.)*.15;

}


/*

	Round Voronoi Borders
	----------------------

    Producing rounded, evenly distributed Voronoi borders via rudimentary alterations to
	the regular Voronoi formula. As an aside, the results are presented in a faux 3D style.
	Essentially, it's a 2D effect.

	Dr2 has been putting up some rounded Voronoi border examples lately, and Abje produced
	a really cool one using a very simple tweak.

	Dr2's variation is fast and nicely distributed,	and as such, translates well to a
	raymarching environment. Abje's tweak can be combined with either IQ or Tomk's line
	distance Voronoi examples to produce really good quality rounded borders - I intend to
	produce an example later that I hope does it justice.

	This example utilizes yet another variation that I put together ages ago. I've outlined
	the method in the Voronoi function - not that it needs much explaining. It does the job
	under then right circumstances and it's reasonably cheap and simple to implement. However,
	for robustness, I'd suggest using one of the aforementioned methods.

	By the way, all variations basically do the same thing, and rely on the idea of
	incorporating a smooth distance metric into a Voronoi-like formula, which IQ wrote about
	in his article on smooth Voronoi.

	I should also mention that Fabrice Neyret incorporated a third order distance to produce
	a rounded border effect, which I used for an example a while back.

	Anyway, just for fun, I like to make 3D looking effects using nothing more than 2D layers.
	In this case, I went for a vector layered kind of aesthetic. For all intents and purposes,
    this example is just a few layers strategically laced together. It's all trickery, so
	there's very little physics involved.

	Basically, I've taken a Voronoi sample, then smoothstepped it in various ways to produce
	the web-like look. I've also taken two extra nearby samples in opposite directions,
	then combined the differences to produce opposing gradients to give highlights, the red
	and blue environmental reflections, etc. There's an offset layer for fake shadowing,
	the function value is used for fake occusion... It's all fake, and pretty simple too. :)


    // Other examples:

    // Faster method, and more evenly distributed.
    Smoothed Voronoi Tunnel - Dr2
	https://www.shadertoy.com/view/4slfWl

	// I like this method, and would like to cover it at some stage.
	Round Voronoi - abje
	https://www.shadertoy.com/view/ldXBDs

	// Smooth Voronoi distance metrics. Not about round borders in particular, but it's
	// the idea from which everything is derived.
	Voronoise - Article, by IQ
	http://iquilezles.org/www/articles/voronoise/voronoise.htm

    // A 3rd order nodal approach - I used it in one of my examples a while back.
	2D trabeculum - FabriceNeyret2
	https://www.shadertoy.com/view/4dKSDV


*/

// Define for a some hacky random cell coloring.
//#define RAND_CELL_COLOR

//#ifdef RAND_CELL_COLOR
vec2 id; // Global unique cell identifier.
float hash(float x) { return fract(sin(x)*43758.5453); } // IQ's float to float hash.
//#endif

// vec2 to vec2 hash.
vec2 hash22(vec2 p) {

    // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked
    // amalgamation I put together, based on a couple of other random algorithms I've
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(41, 289)));
    //return fract(vec2(262144, 32768)*n);

    // Animated.
    p = fract(vec2(262144, 32768)*n);
    // Note the ".333," insted of ".5" that you'd expect to see. When edging, it can open
    // up the cells ever so slightly for a more even spread. In fact, lower numbers work
    // even better, but then the random movement would become too restricted. Zero would
    // give you square cells.
    vec2 anim1 = sin( p*6.2831853 + morph_time)*.333 + .333;
    vec2 anim2 = sin( p*6.2831853 + morph_time*2.)*(cos( p*6.2831853 + morph_time*.5)*.3 + .5)*.45 + .5;

    return mix(anim1, anim2, morph_type);
}

// IQ's smooth minimum function.
float smin(float a, float b, float k){

    float h = clamp(.5 + .5*(b - a)/k, 0., 1.);
    return mix(b, a, h) - k*h*(1. - h);
}

// Commutative smooth minimum function. Provided by Tomkh and taken from
// Alex Evans's (aka Statix) talk:
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin2(float a, float b, float r)
{
   float f = max(0., 1. - abs(b - a)/r);
   return min(a, b) - r*.25*f*f;
}

// IQ's exponential-based smooth minimum function. Unlike the polynomial-based
// smooth minimum, this one is associative and commutative.
float sminExp(float a, float b, float k)
{
    float res = exp(-k*a) + exp(-k*b);
    return -log(res)/k;
}




// 2D 2nd-order Voronoi: Obviously, this is just a rehash of IQ's original. I've tidied
// up those if-statements. Since there's less writing, it should go faster. That's how
// it works, right? :)
//
// This is exactly like a regular Voronoi function, with the exception of the smooth
// distance metrics.
float Voronoi(in vec2 p){

    // Partitioning the grid into unit squares and determining the fractional position.
	vec2 g = floor(p), o; p -= g;

    // "d.x" and "d.y" represent the closest and second closest distances
    // respectively, and "d.z" holds the distance comparison value.
	vec3 d = vec3(8.); // 8., 2, 1.4, etc.



    // A 4x4 grid sample is required for the smooth minimum version.
	for(int j = -1; j <= 2; j++){
		for(int i = -1; i <= 2; i++){

			o = vec2(i, j); // Grid reference.
             // Note the offset distance restriction in the hash function.
            o += hash22(g + o) - p; // Current position to offset point vector.

            // Distance metric. Unfortunately, the Euclidean distance needs
            // to be used for clean equidistant-looking cell border lines.
            // Having said that, there might be a way around it, but this isn't
            // a GPU intensive example, so I'm sure it'll be fine.
			d.z = length(o);

            // Hacked in random ID. There'd be smarter ways to do this.
            //#ifdef RAND_CELL_COLOR
            if(d.z<d.x) id = g + vec2(i, j);
            //#endif

            // Up until this point, it's been a regular Voronoi example. The only
            // difference here is the the mild smooth minimum's to round things
            // off a bit. Replace with regular mimimum functions and it goes back
            // to a regular second order Voronoi example.
            d.y = max(d.x, smin(d.y, d.z, .4)); // Second closest point with smoothing factor.

            d.x = smin(d.x, d.z, .2); // Closest point with smoothing factor.

            // Based on IQ's suggestion - A commutative exponential-based smooth minimum.
            // This algorithm is just an approximation, so it doesn't make much of a difference,
            // but it's here anyway.
            //d.y = max(d.x, sminExp(d.y, d.z, 10.)); // Second closest point with smoothing factor.
            //d.x = sminExp(d.x, d.z, 20.); // Closest point with smoothing factor.


		}
	}

    // Return the regular second closest minus closest (F2 - F1) distance.
    return d.y - d.x;

}

// MS: Added this, pulled from github: https://gist.github.com/mairod/a75e7b44f68110e1576d77419d608786
vec3 hueShift(vec3 color, float hue) {
    const vec3 k = vec3(0.57735, 0.57735, 0.57735);
    float cosAngle = cos(hue);
    return vec3(color * cosAngle + cross(k, color) * sin(hue) + k * dot(k, color) * (1.0 - cosAngle));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Screen coordinates.
	vec2 uv = (fragCoord.xy - RENDERSIZE.xy*.5)/ RENDERSIZE.y;

    if (_exists(syn_UserImage)) uv = mix(uv, _loadUserImage().xy, image_mix);
    //
    // Mild, simplistic fisheye effect.
    uv = mix(uv, uv*(.9 + length(uv)*.2), fish_eye);
    //
    // Scrolling action.
    uv -= fly_time*vec2(1, .25)/8.;

    // The function samples: Six 4x4 grid Voronoi function calls. That amount of work would
    // break an old computer, but it's nothing for any reasonably modern GPU.
    //
    // Base function value.
    float c = Voronoi(uv*alt_voronoi);
    //#ifdef RAND_CELL_COLOR
    vec2 svID = id; // Save the global cell ID.
    //#endif
    // Nearby samples to the bottom right and bottom left.
    float c2 = Voronoi(uv*5. - .002);
    float c3 = Voronoi(uv*5. + .002);
    // A more distant sample - used to fake a shadow and highlight.
    float c4 = Voronoi(uv*5. + vec2(.7, 1)*.2);
    // Slight warped finer detailed (higher frequency) samples.
    float c15 = Voronoi(uv*15. + c);
    float c45 = Voronoi(uv*45. + c*2.);


    // Shading the Voronoi pattern.
    //
    // Base shading and a mild spotty pattern.
    vec3 col = vec3(c*c)*(.9 + (c15 - smoothstep(.2, .3, c15))*(.2+spots));

    //
    // Mixing in some finer cloudy detail.
    float sv = c15*.66 + (1. - c45)*(.34+overlay); // Finer overlay pattern.
    col = col*.8 + sv*sqrt(sv)*.4; // Mix in a little of the overlay.


    // Simple pixelated grid overlay for a mild pixelated effect.
    vec2 sl = mod(fragCoord, 2.);
    //
    // It looks more complicated than it is. Mildly darken every second vertical
    // and horizontal pixel, and mildly lighten the others.
    col *= 4.*(1. + step(1., sl.x))*(1. + step(1., sl.y))/9.;


    //#ifdef RAND_CELL_COLOR
    // Random cell coloring.
    if(hash(svID.x*113. + svID.y)>.5) col *= vec3(1.4, .85, 1.25);
    //#endif


    // Adding a red highlight to the bottom left and a blue highlight to the top right. The
    // end result is raised bubbles with environmental reflections. All fake, of course...
    // Having said that, there is a little directional derivative science behind it.
    float b1 = max((c2 - c)/.002, 0.); // Gradient (or bump) factor 1.
    float b2 = max((c3 - c)/.002, 0.); // Gradient (or bump) factor 2.
    //
    // A touch of deep red and blue, with a bit of extra specularity.
    col += vec3(1, .0, .0)*b2*b2*b2*.45 + vec3(0, .0, 1)*b1*b1*b1*.7;
    //
    // Slightly more mild orange and torquoise with less specularity.
    col += vec3(1, .6, .4)*b2*b2*.3 + vec3(.4, .2, 1)*b1*b1*.3;

    // Distant sampled overlay for a shadowy highlight effect. Made up. There'd be better ways.
    float bord2 = smoothstep(0., fwidth(c4)*3., c4 - .1);
    col = max(col + (1.-bord2)*.25, 0.);

    // The web-like overlay. Tweaked to look a certain way.
    float bord3 = smoothstep(0., fwidth(c)*3., c - .1) - smoothstep(0., fwidth(c)*2., c - .08);
    col *= 1. + bord3*1.5;

    // Another darker patch overlay to give a shadowy reflected look. Also made up.
    float sh = max(c4 - c, 0.);
	col *= (1. - smoothstep(0.015, .05, sh)*.4);

    // For some reason, I wanted a bit more shadow down here... I'm sure I had my reasons. :)
    col -= (1.-bord2)*.1;



    // Smoothstepping the original function value, then multiplying for oldschool, fake occlusion.
    col *= smoothstep(0., .15, c)*.85 + .15 + syn_MidHighHits*.1*reactive_lights;

    // Postprocessing. Mixing in bit of ramped up color value to bring the color out more.
    col = mix(col, pow(max(col, 0.), vec3(4)), .333);

    col -= syn_BassLevel*reactive_lights*.2;

    // Rought gamma correction and screen presentation.
	fragColor = vec4(sqrt(col), 1);

	return fragColor;
 }



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
