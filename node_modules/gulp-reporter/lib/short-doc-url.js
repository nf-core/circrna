'use strict';
const inGFW = require('in-gfw');
const axios = require('axios');
const ci = require('ci-info');
const locale = require('./locale');

const isInGFW = ci.isCI && ci.name ? Promise.resolve(locale === 'zh_CN') : inGFW('goo.gl', 't.cn');

const shorturlCache = isInGFW.then(inGFW => {
	return require(inGFW ? './shorturl_cn.json' : './shorturl.json');
});

function shortDocUrl (error) {
	if (!error.doc) {
		return error;
	}
	return shortUrl(error.doc).catch(ex => {
		//
	}).then(shortUrl => {
		if (shortUrl) {
			error.docShort = shortUrl;
		}
		return error;
	});
}

function shortUrl (url) {
	return shorturlCache.then(shorturlCache => {
		const urlLowerCase = url.toLowerCase();
		if (shorturlCache[urlLowerCase]) {
			return shorturlCache[urlLowerCase];
		}
		return isInGFW.then(inGFW => (
			inGFW
				? axios.get(
					`https://api.t.sina.com.cn/short_url/shorten.json?source=3271760578&url_long=${encodeURIComponent(url)}`
				).then(result => (
					result.data[0].url_short.replace(/^https?/i, 'https')
				))
				: axios.post('https://www.googleapis.com/urlshortener/v1/url?key=AIzaSyACqNSi3cybDvDfWMaPyXZEzQ6IeaPehLE', {
					'longUrl': url,
				}).then(result =>
					result.data.id
				)
		));
	});
}

module.exports = function (errors) {
	return Promise.all(errors.map(shortDocUrl));
};
