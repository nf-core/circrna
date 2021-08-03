'use strict';

const BigNumber = require('bignumber.js');

const TimeFrame = require('./TimeFrame');
const constants = require('./constants');

const units = {};
units.y = units.year = units.years = constants.SECONDS_IN_YEAR;
units.w = units.week = units.weeks = constants.SECONDS_IN_WEEK;
units.d = units.day = units.days = constants.SECONDS_IN_DAY;
units.h = units.hour = units.hours = constants.SECONDS_IN_HOUR;
units.m = units.minute = units.minutes = constants.SECONDS_IN_MINUTE;
units.s = units.second = units.seconds = constants.ONE_SECOND;
units.ms = units.millisecond = units.milliseconds = constants.SECONDS_IN_MILLISECOND;

const toTime = (text) => {
    /**
     * Matching the tested string from start to end - /^ .. $/
     * Required - starting with at least one digit - \d+
     * Optional - decimal value - (.[\d]+)?
     * Optional - one or more whitespaces between the number and unit - (\s+)?
     * Required - letters representing time unit (e.g.: years, day, h, s) - [a-zA-Z]+ 
     */
    const matchPattern = /^\d+(.[\d]+)?(\s+)?[a-zA-Z]+$/;

    const splitted = text.toLowerCase()
        //Replace double whitespaces with single whitespace
        .replace(/\s+/gi, ' ')
        //Remove all commas
        .replace(/[,]/g, '')
        //Used for supporting strings such as '30 Years, 2 Hours'
        .replace(/([\d])\s([\w])/gi, '$1$2')
        //Split array to number-type combinations
        .split(' ');

    const errors = splitted.some(v => {
        return v.match(matchPattern) === null;
    });

    if (errors) {
        throw new Error('Invalid format');
    }

    const seconds = splitted
        .filter(v => {
            return v.match(matchPattern);
        })
        .map(v => {
            const numeric = v.match(/\d+(.[\d]+)?(\s+)?/)[0];
            return [parseFloat(v.substr(0, numeric.length)), v.substr(numeric.length).toLowerCase()];
        })
        .map(v => {
            if ( !units.hasOwnProperty(v[1]) ) {
                throw new Error(`Invalid unit: ${v[1]}`);
            }
            return new BigNumber(v[0]).times(units[v[1]]);
        })
        .reduce((prev, curr) => {
            return prev.plus(curr);
        }, new BigNumber(0));

    return new TimeFrame(seconds);
};

//Constructors
toTime.fromMilliseconds = (milliseconds) => {
    return new TimeFrame(milliseconds * constants.SECONDS_IN_MILLISECOND);
};

toTime.fromSeconds = (seconds) => {
    return new TimeFrame(seconds);
};

toTime.fromMinutes = (minutes) => {
    return new TimeFrame(minutes * constants.SECONDS_IN_MINUTE);
};

toTime.fromHours = (hours) => {
    return new TimeFrame(hours * constants.SECONDS_IN_HOUR);
};

toTime.fromDays = (days) => {
    return new TimeFrame(days * constants.SECONDS_IN_DAY);
};

toTime.fromWeeks = (weeks) => {
    return new TimeFrame(weeks * constants.SECONDS_IN_WEEK);
};

toTime.fromYears = (years) => {
    return new TimeFrame(years * constants.SECONDS_IN_YEAR);
};

module.exports = toTime;