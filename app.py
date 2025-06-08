from flask import Flask, render_template, request, jsonify
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
from astropy import units as u
from datetime import datetime, timedelta
import pytz
import numpy as np
import requests

app = Flask(__name__)

CONSTELLATIONS = {
    "corona_borealis": {
        "name": "Corona Borealis",
        "star": "Alphecca (α Coronae Borealis)",
        "ra": "15h34m41s",
        "dec": "+26d42m53s",
    },
    "camelopardalis": {
        "name": "Camelopardalis",
        "star": "β Camelopardalis",
        "ra": "05h03m25s",
        "dec": "+60d26m32s",
    },
    "coma_berenices": {
        "name": "Coma Berenices",
        "star": "β Comae Berenices",
        "ra": "13h11m52s",
        "dec": "+27d52m41s",
    },
}


def geocode_address(address):
    """Convert address to latitude/longitude using OpenStreetMap Nominatim API"""
    try:
        base_url = "https://nominatim.openstreetmap.org/search"
        params = {"q": address, "format": "json", "limit": 1, "addressdetails": 1}
        headers = {"User-Agent": "ZenithCalculator/1.0"}

        response = requests.get(base_url, params=params, headers=headers, timeout=5)
        response.raise_for_status()

        data = response.json()
        if not data:
            return {"success": False, "error": "Address not found"}

        result = data[0]
        return {
            "success": True,
            "latitude": float(result["lat"]),
            "longitude": float(result["lon"]),
            "display_name": result["display_name"],
        }

    except requests.RequestException as e:
        return {"success": False, "error": "Geocoding service unavailable"}
    except (KeyError, ValueError, IndexError) as e:
        return {"success": False, "error": "Invalid response from geocoding service"}


def calculate_zenith_time(
    date, latitude, longitude, constellation_key="corona_borealis"
):
    """Calculate when constellation's brightest star reaches maximum altitude for given location and date"""
    try:
        # Get constellation data
        if constellation_key not in CONSTELLATIONS:
            raise ValueError(f"Unknown constellation: {constellation_key}")

        constellation = CONSTELLATIONS[constellation_key]

        # Create coordinates for the constellation's star
        star = SkyCoord(ra=constellation["ra"], dec=constellation["dec"], frame="icrs")

        # Create observer location
        location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg)

        # Create time array for the given date (24 hours)
        utc_date = datetime.strptime(date, "%Y-%m-%d")
        times = Time(
            [utc_date + timedelta(hours=i / 4) for i in range(96)]
        )  # 15-minute intervals

        # Calculate altitude and azimuth
        altaz_frame = AltAz(obstime=times, location=location)
        star_altaz = star.transform_to(altaz_frame)

        # Find maximum altitude
        max_alt_idx = np.argmax(star_altaz.alt)
        max_altitude = star_altaz.alt[max_alt_idx].deg
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
            "constellation_name": constellation["name"],
            "star_name": constellation["star"],
        }

    except Exception as e:
        raise e
        return {"success": False, "error": str(e)}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/geocode", methods=["POST"])
def geocode():
    try:
        address = request.form.get("address", "").strip()
        if not address:
            return jsonify({"success": False, "error": "Address is required"})

        result = geocode_address(address)
        return jsonify(result)

    except Exception as e:
        return jsonify({"success": False, "error": "Geocoding failed"})


@app.route("/calculate", methods=["POST"])
def calculate():
    try:
        date = request.form["date"]
        latitude = float(request.form["latitude"])
        longitude = float(request.form["longitude"])
        constellation = request.form.get("constellation", "corona_borealis")

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

        result = calculate_zenith_time(date, latitude, longitude, constellation)
        return jsonify(result)

    except ValueError as e:
        return jsonify({"success": False, "error": "Invalid input values"})
    except Exception as e:
        raise e
        return jsonify({"success": False, "error": str(e)})


if __name__ == "__main__":
    app.run(debug=True)
