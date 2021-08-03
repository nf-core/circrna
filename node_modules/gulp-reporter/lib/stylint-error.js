'use strict';
const LintError = require('./lint-error');
const STYLLINT_SEVERITY_MAP = {
	'warning': 'warn',
};

/**
 * stylint错误对象
 *
 * @class StyLintError
 * @extends {LintError}
 */
class StyLintError extends LintError {
	/**
	 * Creates an instance of StyLintError.
	 *
	 * @param {Object} message stylint的原始error对象
	 *
	 * @memberOf StyLintError
	 */
	constructor (message) {
		const error = {
			severity: 'error',

			// 报错插件
			plugin: 'StyLint',
		};
		message.replace(/^(.+?):\s*(.+?)$/gm, (s, key, value) => {
			if (/^(?:error|warn(?:ing)?|info)$/i.test(key)) {
				key = key.toLowerCase();
				error.severity = STYLLINT_SEVERITY_MAP[key] || key;
				key = 'message';
			} else if (/^file(?:Name)?$/i.test(key)) {
				key = 'fileName';
			} else if (/^line(?:Number)?$/i.test(key)) {
				key = 'lineNumber';
				value = +value.replace(/^(\d+)(?::\s*)?(.*?)$/, '$1');
				if (RegExp.$2) {
					error.source = RegExp.$2;
				}
			} else {
				key = key.toLocaleLowerCase();
			}
			error[key] = value;
		});

		super(error);
	}
}

module.exports = StyLintError;
