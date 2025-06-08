from flask import Flask, render_template, request, jsonify
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
from astropy import units as u
from datetime import datetime, timedelta
import pytz
import numpy as np

app = Flask(__name__)

ALPHECCA_RA = "15h34m41s"
ALPHECCA_DEC = "+26d42m53s"


def calculate_zenith_time(date, latitude, longitude):
    """Calculate when Alphecca reaches maximum altitude for given location and date"""
    try:
        # Create coordinates for Alphecca
        alphecca = SkyCoord(ra=ALPHECCA_RA, dec=ALPHECCA_DEC, frame="icrs")

        # Create observer location
        location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg)

        # Create time array for the given date (24 hours)
        utc_date = datetime.strptime(date, "%Y-%m-%d")
        times = Time(
            [utc_date + timedelta(hours=i / 4) for i in range(96)]
        )  # 15-minute intervals

        # Calculate altitude and azimuth
        altaz_frame = AltAz(obstime=times, location=location)
        alphecca_altaz = alphecca.transform_to(altaz_frame)

        # Find maximum altitude
        max_alt_idx = np.argmax(alphecca_altaz.alt)
        max_altitude = alphecca_altaz.alt[max_alt_idx].deg
        zenith_time_utc = times[max_alt_idx].datetime

        # Convert to local timezone (approximate based on longitude)
        tz_offset = int(longitude / 15)  # Rough timezone calculation
        local_time = zenith_time_utc + timedelta(hours=tz_offset)

        # Check if observable (not during daylight)
        sun_altaz = get_sun(times[max_alt_idx]).transform_to(altaz_frame[max_alt_idx])
        sun_altitude = sun_altaz.alt.deg
        is_observable = sun_altitude < -18  # Astronomical twilight

        # Check if it actually reaches zenith (>80 degrees is considered "near zenith")
        reaches_zenith = max_altitude > 80

        return {
            "success": True,
            "zenith_time": local_time.strftime("%H:%M:%S"),
            "max_altitude": round(max_altitude, 1),
            "is_observable": bool(is_observable),
            "reaches_zenith": bool(reaches_zenith),
            "sun_altitude": round(sun_altitude, 1),
        }

    except Exception as e:
        raise e
        return {"success": False, "error": str(e)}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/calculate", methods=["POST"])
def calculate():
    try:
        date = request.form["date"]
        latitude = float(request.form["latitude"])
        longitude = float(request.form["longitude"])

        # Validate inputs
        if not (-90 <= latitude <= 90):
            return jsonify(
                {
                    "success": False,
                    "error": "Latitude must be between -90 and 90 degrees",
                }
            )
        if not (-180 <= longitude <= 180):
            return jsonify(
                {
                    "success": False,
                    "error": "Longitude must be between -180 and 180 degrees",
                }
            )

        result = calculate_zenith_time(date, latitude, longitude)
        return jsonify(result)

    except ValueError as e:
        return jsonify({"success": False, "error": "Invalid input values"})
    except Exception as e:
        raise e
        return jsonify({"success": False, "error": str(e)})


if __name__ == "__main__":
    app.run(debug=True)
