import ee

def add_band(image):
    """
    Add total precipitable water values and index for the LUT of SMW algorithm coefficients to the image.

    Parameters:
    - image (gee.Image): Image for which to interpolate the TPW data. Needs the 'system:time_start' image property.

    Returns:
    - gee.Image: Image with added 'TPW' and 'TPWpos' bands.
    """

    date = image.get('system:time_start')
    year, month, day = date.year(), date.month(), date.day()
    date1 = gee.Date.from_ymd(year, month, day)
    date2 = date1.advance(1, 'days')

    def datedist(img):
        return img.set('DateDist', abs(img.get('system:time_start') - date.millis()))

    TPWcollection = (gee.ImageCollection('NCEP_RE/surface_wv')
                     .filter_date(date1.format('yyyy-MM-dd'), date2.format('yyyy-MM-dd'))
                     .map(datedist))

    closest = TPWcollection.sort('DateDist').to_list(2)

    tpw1 = closest.get(0).select('pr_wtr') if closest.size() else gee.Image.constant(-999.0)
    tpw2 = closest.get(1).select('pr_wtr') if closest.size() > 1 else tpw1

    time1 = tpw1.get('DateDist') / 21600000 if closest.size() else 1.0
    time2 = tpw2.get('DateDist') / 21600000 if closest.size() > 1 else 0.0

    tpw = tpw1.expression('tpw1*time2+tpw2*time1', {
        'tpw1': tpw1,
        'time1': time1,
        'tpw2': tpw2,
        'time2': time2
    }).clip(image.geometry())

    pos = tpw.expression(
        "value = (TPW>0 && TPW<=6) ? 0" +
        ": (TPW>6 && TPW<=12) ? 1" +
        ": (TPW>12 && TPW<=18) ? 2" +
        ": (TPW>18 && TPW<=24) ? 3" +
        ": (TPW>24 && TPW<=30) ? 4" +
        ": (TPW>30 && TPW<=36) ? 5" +
        ": (TPW>36 && TPW<=42) ? 6" +
        ": (TPW>42 && TPW<=48) ? 7" +
        ": (TPW>48 && TPW<=54) ? 8" +
        ": (TPW>54) ? 9" +
        ": 0", {'TPW': tpw}).clip(image.geometry())

    withTPW = image.add_bands(tpw.rename('TPW')).add_bands(pos.rename('TPWpos'))

    return withTPW
