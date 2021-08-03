'use strict';
module.exports = (() => {
	try {
		/* global navigator */
		return navigator.language.replace(/-/g, '_');
	} catch (ex) {
		const osLocale = require('os-locale');
		return osLocale.sync();
	}
})();
