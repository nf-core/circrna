'use strict';
const locale = require('./locale');

const I18N_CACHE = {

};

function i18n (error) {
	if (locale === 'en_US') {
		return error;
	}
	const key = error.plugin.toLowerCase();
	if (!(key in I18N_CACHE)) {
		try {
			I18N_CACHE[key] = require(`./${key}_${locale}`);
		} catch (ex) {
			I18N_CACHE[key] = false;
		}
	}
	if (I18N_CACHE[key]) {
		return I18N_CACHE[key](error);
	}
}

/**
 * 统一的，错误对象
 *
 * @class LintError
 * @extends {Error}
 */
class LintError extends Error {
	/**
	 * Creates an instance of LintError.
	 *
	 *
	 * @memberOf LintError
	 */
	constructor () {
		super();
		const args = Array.from(arguments);
		args.unshift(this);
		Object.assign.apply(Object, args);
		i18n(this);
		this.name = this.constructor.name;
		delete this.stack;
	}

	/**
	 * console.log会自动调用的函数
	 *
	 * @returns {String} 模拟Error对象的格式的错误信息的字符串
	 *
	 * @memberOf LintError
	 */
	inspect () {
		const message = `${this.name}: ${this.message} (${[
			this.plugin,
			this.rule,
		].filter(Boolean).join(' ')})`;
		let stack = [
			[
				this.fileName,
				this.lineNumber,
				this.columnNumber,
			].filter(Boolean).join(':'),
		];
		if (this.doc) {
			stack.push(this.doc);
		}
		stack = stack.map((uri) => {
			return '    at ' + uri;
		});
		stack.unshift(message);
		return stack.join('\n');
	}
	get stack () {
		return this.inspect();
	}
}

module.exports = LintError;
