"use strict";

let result = null;
function async (mod, parseValue) {
	let lazyResult;

	const promise = () => (
		new Promise((resolve, reject) => {
			process.nextTick(() => {
				if (result != null) {
					return resolve(result);
				}
				mod.async((error, value) => {
					if (error) {
						reject(error);
					} else {
						resolve(parseValue(value));
					}
				});
			});
		})
	);

	function lazyPromise () {
		if (!lazyResult) {
			lazyResult = promise();
		}
		return lazyResult;
	}

	return function (callback) {
		if (callback) {
			lazyPromise().then(() => {
				callback(null, result);
			}, callback);
		} else {
			return lazyPromise();
		}
	};
}

function sync (mod, parseValue) {
	return function () {
		if (result != null) {
			return result;
		}
		let resultError, resultValue;
		mod.sync((error, value) => {
			resultError = error;
			resultValue = value;
		});
		if (resultError) {
			throw resultError;
		} else {
			return parseValue(resultValue);
		}
	};
}

function mod (name, pattern) {
	function parseValue (value) {
		result = pattern.test(value);
		return result;
	}

	const mod = require("./" + name);
	module.exports = {
		async: async(mod, parseValue),
		sync: sync(mod, parseValue),
	};
}

if (process.platform === "win32" || require("is-wsl")) {
	mod("geo", /^0*45$/);
} else {
	mod("timezone", /(^|\/)(?:PRC|Asia\/(?:Shanghai|Chongqing|Urumqi|Beijing))$/i);
}
