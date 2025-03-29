#include "colormaps.c"

double color_amplitude(double value, double scale, int time)
/* transforms the wave amplitude into a double in [-1,1] to feed into color scheme */
{
    return(tanh(SLOPE*value/scale)*exp(-((double)time*ATTENUATION)));
}

double color_amplitude_asym(double value, double scale, int time)
/* transforms the wave amplitude into a double in [-1,1] to feed into color scheme */
{
    return(2.0*tanh(SLOPE*value/scale)*exp(-((double)time*ATTENUATION)) - 1.0);
}

double color_amplitude_linear(double value, double scale, int time)
/* transforms the wave amplitude into a double in [-1,1] to feed into color scheme */
{
    return(SLOPE*value/scale);
}

void color_scheme(int scheme, double value, double scale, int time, double rgb[3]) /* color scheme */
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case (C_LUM):
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case (C_HUE):
        {
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = 0.5;
            hue = HUEMEAN + amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case (C_ONEDIM):
        {
            amplitude = color_amplitude(value, scale, time);
            amp_to_rgb(0.5*(1.0 + amplitude), rgb);
            break;
        }
        case (C_ONEDIM_LINEAR):
        {
            amplitude = color_amplitude_linear(value, scale, time);
            if (amplitude > 1.0) amplitude -= 1.0;
            else if (amplitude < 0.0) amplitude += 1.0;
            amp_to_rgb(amplitude, rgb);
            break;
        }
        case (C_BASIC_LINEAR):
        {
            amplitude = color_amplitude_linear(value, scale, time);
            if (amplitude > 1.0) amplitude -= 1.0;
            else if (amplitude < 0.0) amplitude += 1.0;
            amp_to_rgb(amplitude, rgb);
            break;
        }
    }
}


void color_scheme_lum(int scheme, double value, double scale, int time, double lum, double rgb[3]) /* color scheme */
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case (C_LUM):
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case (C_HUE):
        {
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = lum;
            hue = HUEMEAN + amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
    }
}

void color_scheme_asym(int scheme, double value, double scale, int time, double rgb[3]) /* color scheme */
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case (C_LUM):
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case (C_HUE):
        {
            r = 0.9;
            amplitude = color_amplitude_asym(value, scale, time);
            y = 0.5;
            hue = HUEMEAN + 0.8*amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case (C_ONEDIM):
        {
            amplitude = color_amplitude(value, scale, time);
            amp_to_rgb(amplitude, rgb);
            break;
        }
    }
}

void color_scheme_palette(int scheme, int palette, double value, double scale, int time, double rgb[3]) /* color scheme */
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case (C_LUM):
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb_palette(hue, r, y, rgb, palette);
            break;
        }
        case (C_HUE):
        {
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = 0.5;
            hue = HUEMEAN + amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb_palette(hue, r, y, rgb, palette);
            break;
        }
        case (C_ONEDIM):
        {
            amplitude = color_amplitude(value, scale, time);
            amp_to_rgb_palette(0.5*(1.0 + amplitude), rgb, palette);
            break;
        }
        case (C_ONEDIM_LINEAR):
        {
            amplitude = color_amplitude_linear(value, scale, time);
            if (amplitude > 1.0) amplitude -= 1.0;
            else if (amplitude < 0.0) amplitude += 1.0;
            amp_to_rgb_palette(amplitude, rgb, palette);
            break;
        }
        case (C_BASIC_LINEAR):
        {
            if (value > 1.0) value -= 1.0;
            else if (value < 0.0) value += 1.0;
            amp_to_rgb_palette(value, rgb, palette);
            break;
//             if (value > 2.0) value -= 2.0;
//             else if (value < 0.0) value += 2.0;
//             amp_to_rgb_palette(value, rgb, palette);
//             break;
        }
    }
}

void color_scheme_asym_palette(int scheme, int palette, double value, double scale, int time, double rgb[3]) /* color scheme */
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case (C_LUM):
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb_palette(hue, r, y, rgb, palette);
            break;
        }
        case (C_HUE):
        {
            r = 0.9;
            amplitude = color_amplitude_asym(value, scale, time);
            y = 0.5;
            hue = HUEMEAN + 0.8*amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb_palette(hue, r, y, rgb, palette);
            break;
        }
        case (C_ONEDIM):
        {
            amplitude = color_amplitude(value, scale, time);
            amp_to_rgb_palette(amplitude, rgb, palette);
            break;
        }
    }
}


