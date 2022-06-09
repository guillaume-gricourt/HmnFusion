from typing import Any, Dict

import jinja2
from hmnfusion import _version


class Report(object):
    """Report builds reports for HmnFusion"""

    @classmethod
    def select_template(cls, template: str) -> jinja2.Template:
        """Select template in the package

        Parameters
        ----------
        template: str
            Path of the template to select

        Return
        ------
        jinja2.Template
            the selected template
        """
        env = jinja2.Environment(
            loader=jinja2.PackageLoader(_version.__app_name__),
            autoescape=jinja2.select_autoescape(),
        )
        return env.get_template(template)

    @classmethod
    def to_html(
        cls, path: str, template: jinja2.Template, context: Dict[Any, Any]
    ) -> None:
        """Write context and template into a file.

        Parameters
        ----------
        path: str
            Path of an output file (xlsx, excel format)
        template: jinja2.Template
            Template into write
        context: dict
            Data to transfert into the template

        Return
        ------
        None
        """
        msg = template.render(context=context)
        with open(path, "w") as fod:
            fod.write(msg)
