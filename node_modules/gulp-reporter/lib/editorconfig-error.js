'use strict';
const LintError = require('./lint-error');
/**
 * EditorConfig错误对象
 *
 * @class EditorConfigError
 * @extends {LintError}
 */
class EditorConfigError extends LintError {
	/**
	 * Creates an instance of EditorConfigError.
	 *
	 * @param {Object} error eclint原生错误对象
	 *
	 * @memberOf EditorConfigError
	 */
	constructor (error) {
		super({
			// EditorConfig无警告，错误等级全部算错误
			severity: 'error',

			// 报错插件
			plugin: 'EditorConfig',

			// 文档
			doc: error.rule && `https://github.com/editorconfig/editorconfig/wiki/EditorConfig-Properties#${error.rule}`,
		}, error);
	}
}

module.exports = EditorConfigError;
